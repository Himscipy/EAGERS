function Locked = CheckRampRates(QP,Locked,OptimalState,dt)
global Plant
%Premise, it may not be possible to follow the output from the step-by-step
%dispatch due to ramp rate limitations. The ramp rates will limit the 
%cumulative energy generation. Storage or spare capacity in other generators 
%can make up the difference without turning something on. 
%Otherwise, something must turn on.
nG = length(Locked(1,:));
networkNames = fieldnames(QP.Organize.Balance);
%% make sure it can shut down in time from initial condition
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1
        Locked(OptimalState(:,i)==0,i)=false;%Identify when optimal dispatch has things off
        if OptimalState(1,i)>0 && ~all(Locked(:,i))
            r = QP.Organize.Ramping(i)+1;
            D = OptimalState(1,i);
            t = 1;
            while D>0
                D = D - QP.b(r);
                if D>0 && ~Locked(t+1,i)
                    Locked(t+1,i) = true;
                end
                t = t+1;
                r = r+QP.Organize.t1ineq;
            end
        end
    end
end

%% Make sure the loss of energy due to ramping constraints at starts & stops, can be made up by other generators or storage
for net = 1:1:length(networkNames)
    inc = false(nG,1);
    stor = false(nG,1);
    util = [];
    out = Plant.subNet.(networkNames{net}).abbreviation;
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'Solar') 
            %%don't include
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            if isfield(Plant.Generator(i).QPform.output,out) || (strcmp(out,'DC') && isfield(Plant.Generator(i).QPform.output,'E'))
                util = i;
            end
        elseif ismember(Plant.Generator(i).Type,{'Thermal Storage';'Electric Storage';}) 
            if isfield(Plant.Generator(i).QPform.output,out)
                stor(i) = true;
            end
        elseif strcmp(Plant.Generator(i).Type,'Chiller') 
            if strcmp(out,'C')
                inc(i) = true;
            end
        elseif isfield(Plant.Generator(i).QPform.output,out)
            inc(i) = true;
        end
    end
    if any(stor) %&& any(StoredEnergy>0)
        Locked = checkStorageCapacity(QP,OptimalState,stor,Locked,inc,out,dt);
    else
        Locked = turnSomethingOn(QP,OptimalState,Locked,out,networkNames{net},inc,util,dt);
    end
end
end%ends function CheckRampRates

function Locked = checkStorageCapacity(QP,OptimalState,stor,Locked,inc,out,dt)
%start something earlier, or leave something on, if there is not enough
%storage capacity to makeup the difference between the optimal dispatch and
%the maximum output given when the generator turns on/off
global Plant
nG = length(inc);
nS = length(dt);
StoredEnergy = zeros(nS,1);
ChargeEff = zeros(nG,1);
for i = 1:1:nG
    if stor(i)%find the cumulative stored energy of this type
        buff = Plant.Generator(i).QPform.Stor.UsableSize*(Plant.Generator(i).VariableStruct.Buffer/100);
        StoredEnergy = StoredEnergy + max(0,(OptimalState(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff);
        ChargeEff(i) = Plant.Generator(i).QPform.Stor.ChargeEff;
    end
end
avgChargeEff = mean(nonzeros(ChargeEff));
[MaxOut,Constraint] = GenLimit(OptimalState,Locked,dt);
if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
    [Disp,HeatConsumed] = heatDisp(OptimalState,Locked,QP,dt);
else
    Disp = OptimalState;
end
for i = 1:1:nG
    if inc(i)
        starts = nonzeros((1:nS)'.*((Locked(2:end,i)-Locked(1:nS,i))>0));
        k = 1;
        while ~isempty(starts) && k<=length(starts)
            t = starts(k);
            D = zeros(nS,1);
            S = zeros(nS,1);
            newOn = [];
            while t<=nS && MaxOut.(out)(t+1,i)<Disp(t+1,i)
                if t>starts(k)
                    D(t) = D(t-1);
                end
                D(t) = D(t)+Disp(t+1,i)-MaxOut.(out)(t+1,i);%cumulative difference between desired output and feasible output
                t = t+1;
            end
            t_end = max(1,t-1);
            for t = 1:1:t_end%cumulative spare capacity in other generators
                if t>1
                    S(t) = S(t-1);
                end
                if strcmp(out,'H')
                    for j = 1:1:nG 
                        if inc(j) && j~=i
                            S(t) = S(t) + MaxOut.(out)(t+1,j);
                        end
                        S(t) = S(t) - HeatConsumed(t);
                    end
                else
                    for j = 1:1:nG 
                        if inc(j) %&& j~=i
                            S(t) = S(t) + MaxOut.(out)(t+1,j) - Disp(t+1,j);
                        end
                    end
                end
            end
            if any((D - S) > StoredEnergy) && ~isinf(Constraint(starts(k),i)) && ~(Constraint(starts(k),i)==0) && ~Locked(Constraint(starts(k),i)+1,i)
                %update Locked, MaxOut, and StoredEnergy if the generator was started earlier
                newOn = [Constraint(starts(k),i), newOn];
                Locked(newOn+1,i) = true;
                OldMax = MaxOut.(out)(newOn+1,i);
                [MaxOut,Constraint] = GenLimit(OptimalState,Locked,dt);
                OptimalState(newOn+1,i) = MaxOut.(out)(newOn+1,i);
                if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
                    [Disp,HeatConsumed] = heatDisp(OptimalState,Locked,QP,dt);
                else
                    Disp = OptimalState;
                end
                for j = 1:1:length(newOn)
                    StoredEnergy(newOn(j):end) = StoredEnergy(newOn(j):end) + (MaxOut.(out)(newOn(j)+1,i)-OldMax(j))*avgChargeEff;
                end
                if starts(k)>1 && Locked(starts(k),i) == false %move start and check again
                    starts(k) = starts(k)-1;
                    k = k-1;
                end
            end
            k = k+1;
        end

        stops = nonzeros((1:nS)'.*(Locked(1:nS,i)-(Locked(2:end,i))>0));
        if ~isempty(stops) && stops(1) ==1 %Initial condition ramp-down was taken care of earlier
            if length(stops)>1
                stops = stops(2:end);
            else
                stops = [];
            end
        end
        k = 1;
        while ~isempty(stops) && k<=length(stops)
            D = zeros(nS,1);
            S = zeros(nS,1);
            newOn = [];
            nextOn = starts(starts>stops(k));
            if isempty(nextOn)
                nextOn = nS;
            else
                nextOn = nextOn(1);
            end
            t_startDown = 1;
            while Disp(t_startDown+1,i)<MaxOut.(out)(t_startDown+1,i) && t_startDown<stops(k)
                t_startDown = t_startDown+1;%find first time it exceeds ramping constraint
            end
            for t = 1:1:nextOn %cumulative spare capacity in other generators
                if t>1
                    S(t) = S(t-1);
                    D(t) = D(t-1);
                end
                if t>=t_startDown
                    D(t) = max(0,D(t)+Disp(t+1,i)-MaxOut.(out)(t+1,i));%difference between desired output and feasible output
                end
                
                if strcmp(out,'H')
                    for j = 1:1:nG 
                        if inc(j) && j~=i
                            S(t) = S(t) + MaxOut.(out)(t+1,j);
                        end
                        S(t) = S(t) - HeatConsumed(t);
                    end
                else
                    for j = 1:1:nG 
                        if inc(j) && j~=i
                            S(t) = S(t) + MaxOut.(out)(t+1,j) - Disp(t+1,j);
                        end
                    end
                end
            end
            if any((D - S) > StoredEnergy) && Constraint(stops(k)-1,i)>0 && ~isinf(Constraint(stops(k)-1,i)) && ~Locked(Constraint(stops(k)-1,i)+1,i) 
                newOn = [newOn,Constraint(stops(k)-1,i)];
                Locked(newOn+1,i) = true;
                OldMax = MaxOut.(out)(newOn+1,i);
                [MaxOut,Constraint] = GenLimit(OptimalState,Locked,dt);
                OptimalState(newOn+1,i) = MaxOut.(out)(newOn+1,i);
                if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
                    [Disp,HeatConsumed] = heatDisp(OptimalState,Locked,QP,dt);
                else
                    Disp = OptimalState;
                end
                for j = 1:1:length(newOn)
                    StoredEnergy(newOn(j):end) = StoredEnergy(newOn(j):end) + (MaxOut.(out)(newOn(j)+1,i)-OldMax(j))*avgChargeEff;
                end
                if stops(k)<nS && Locked(stops(k)+1,i) == false %move stop and check again
                    stops(k) = stops(k)+1;
                    k = k-1;
                end
            end
            k = k+1;
        end
    end
end
end%Ends function checkStorageCapacity


function [Locked,feas] = turnSomethingOn(QP,OptimalState,Locked,out,net,inc,util,dt)
global Plant
nG = length(Plant.Generator);
nS = length(dt);
lgen = [];
NetDemand = zeros(nS,1);%% Find net demands
for n = 1:1:length(QP.Organize.Balance.(net)) %run through all the nodes in this network
    req = QP.Organize.Balance.(net)(n);%balance at this node (t = 1)
    req = req:QP.Organize.t1Balances:((nS-1)*QP.Organize.t1Balances+req);%balance at this node (t = 1:nS)
    NetDemand = NetDemand + QP.beq(req);
end
for i = 1:1:nG  %add demands from electric or absorption chillers
    if strcmp(Plant.Generator(i).Type,'Chiller') && ((strcmp(out,'H') && strcmp(Plant.Generator(i).Source,'Heat')) || ((strcmp(out,'E') && strcmp(Plant.Generator(i).Source,'Electricity'))) )
        NetDemand = NetDemand + QP.constDemand.(out).load(:,i).*Locked(2:end,i);
        lgen(end+1) = i;
        n = Plant.Generator(i).QPform.(net).subnetNode;
        req = QP.Organize.Balance.(net)(n);%balance at this node (t = 1)
        req = req:QP.Organize.t1Balances:((nS-1)*QP.Organize.t1Balances+req);%balance at this node (t = 1:nS)
        for t = 1:1:nS
            if Locked(t+1,i)
                C = 0;
                j = 1;
                states = QP.Organize.States{i} + (t-1)*QP.Organize.t1States;
                while j<=length(states) && C<OptimalState(t+1,i)
                    c = min(OptimalState(t+1,i) - C,QP.ub(states(j)));
                    C = C + c;
                    NetDemand(t) = NetDemand(t) - c*QP.Aeq(req(t),states(j));
                    j = j+1;
                end
            end
        end
    end
end
%%No storage to check 
% for i = 1:1:nG  %add demands for charging/discharging storage
%     if any(strcmp(Plant.Generator(i).Type, {'Thermal Storage';'Electric Storage'})) && isfield(Plant.Generator(i).QPform.output,out)
%         for t = 1:1:nS
%             if OptimalState(t,i)<OptimalState(t+1,i)%charging
%                 NetDemand(t) = NetDemand(t) + (OptimalState(t+1,i) - OptimalState(t,i))*dt(t)/Plant.Generator(i).QPform.Stor.ChargeEff;
%             else
%                 NetDemand(t) = NetDemand(t) - (OptimalState(t,i) - OptimalState(t+1,i))*dt(t)*Plant.Generator(i).QPform.Stor.DischEff;
%             end
%         end 
%     end
% end
feas = false;
while ~feas
    [MaxOut,Constraint] = GenLimit(OptimalState,Locked,dt);
    if strcmp(out,'H') %overwrite electric dispatch setting with what the heat dispatch is for CHP generators
        [~,NetDemand] = heatDisp(OptimalState,Locked,QP,dt);
    end
    constrainedGen = sum(MaxOut.(out)(2:end,inc),2);
    if ~isempty(util)
        constrainedGen = constrainedGen + MaxOut.(out)(2:end,util);
    end
    if all(constrainedGen>=(NetDemand-1e-9))%fix due to rounding errors
        feas = true;
    else
        feas = false;
        t = 1;
        while sum(MaxOut.(out)(t+1,inc))>=(NetDemand(t)-1e-9)%fix due to rounding errors
            t = t+1;
        end
        diff = zeros(nG,1);
        if any(Constraint(t,inc)>0 & Constraint(t,inc)<nS+1) %at leat one generator in this category is constrained by a startup or shutdown at this time
            %select apropriate one to start early or keep on later
            for i = 1:1:nG
                if inc(i)
                    if Constraint(t,i)>0 && Constraint(t,i)<nS+1
                        diff(i) = abs(Constraint(t,i)-t);
                    else
                        diff(i) = nS+1;
                    end
                else
                    diff(i) = inf;
                end
            end
            [~,I] = min(diff);
            Locked(Constraint(t,I)+1,I) = true;
        elseif any(~Locked(t+1,inc))
            %turn on something else that has the same type of output
            for i = 1:1:nG
                if inc(i)
                    if ~Locked(t+1,i)
                        if t==1
                            diff(i) = nS-t - max((nS-t:-1:1)'.*Locked(t+2:nS+1,i));
                        elseif t == nS
                            diff(i) = t - max((1:t-1)'.*Locked(2:t,i));
                        else
                            diff(i) = min(t - max((1:t-1)'.*Locked(2:t,i)),nS-t - max((nS-t:-1:1)'.*Locked(t+2:nS+1,i)));
                        end
                    else
                        diff(i) = nS+1;
                    end
                else
                    diff(i) = inf;
                end
            end
            [~,I] = min(diff);
            Locked(t+1,I) = true;
        elseif ~isempty(lgen) && any(Locked(t+1,lgen))
            %turn off something that is a load (e.g. a absorption chiller when balancing the heat
            for i = 1:1:length(lgen)
                if Locked(t+1,lgen(i))
                    Locked(t+1,lgen(i)) = false;
                end
            end
        else
            feas = true;%nothing else it can change by this logic
%             disp('Error in CheckRampRates: cant find a generator to turn on to make feasible')
        end
    end
end
end%Ends function turnSomethingOn

function [Disp,HeatConsumed] = heatDisp(OptimalState,Locked,QP,dt)
global Plant
nG = length(Plant.Generator);
[n,~] = size(OptimalState);
nS = n - 1;
Disp = OptimalState; 
HeatConsumed = zeros(nS,1);
for n = 1:1:length(Plant.subNet.DistrictHeat.nodes) %run through all the nodes in this network
    req = QP.Organize.Balance.DistrictHeat(n);%balance at this node (t = 1)
    req = req:QP.Organize.t1Balances:((nS-1)*QP.Organize.t1Balances+req);%balance at this node (t = 1:nS)
    load = nonzeros(QP.Organize.Demand.DistrictHeat{n}); %loads at this node
    if ~isempty(load) %need this in case there is no field Forecast.Demand
        HeatConsumed = HeatConsumed + QP.beq(req); %multiple demands can be at the same node, or none
    end
end
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type, 'CHP Generator')
        for t = 1:1:nS
            D = 0;
            j = 1;
            if Locked(t+1,i)
                Disp(t+1,i) = -Plant.Generator(i).QPform.constDemand.H;               
            else
                Disp(t+1,i) = 0;
            end
            states = QP.Organize.States{i} + (t-1)*QP.Organize.t1States;
            while j<=length(states)
                d = min(OptimalState(t+1,i) - D,QP.ub(states(j)));
                D = D + d;
                Disp(t+1,i) = Disp(t+1,i) + d*QP.Aeq(req(t),states(j));
                j = j+1;
            end
        end
    elseif strcmp(Plant.Generator(i).Type, 'Heater')

    elseif strcmp(Plant.Generator(i).Type, 'Chiller') && strcmp(Plant.Generator(i).Source,'Heat')
        Disp(:,i) = 0;
        HeatConsumed = HeatConsumed + QP.constDemand.H.load(:,i).*Locked(2:end,i);
        for t = 1:1:nS
            if Locked(t+1,i)
                D = 0;
                j = 1;
                states = QP.Organize.States{i} + (t-1)*QP.Organize.t1States;
                while j<=length(states) && D<OptimalState(t+1,i)
                    d = min(OptimalState(t+1,i) - D,QP.ub(states(j)));
                    D = D + d;
                    HeatConsumed(t) = HeatConsumed(t) - d*QP.Aeq(req(t),states(j));
                    j = j+1;
                end
            end
        end
    elseif strcmp(Plant.Generator(i).Type, 'Thermal Storage')
        if strcmp(Plant.Generator(i).Source,'Heat')
            for t = 1:1:nS
                if OptimalState(t,i)<OptimalState(t+1,i)%charging
                    HeatConsumed(t) = HeatConsumed(t) + (OptimalState(t+1,i) - OptimalState(t,i))*dt(t)/Plant.Generator(i).QPform.Stor.ChargeEff;
                else
                    HeatConsumed(t) = HeatConsumed(t) - (OptimalState(t,i) - OptimalState(t+1,i))*dt(t)*Plant.Generator(i).QPform.Stor.DischEff;
                end
            end  
        else
            Disp(:,i) = 0;
        end
    else
        Disp(:,i) = 0;
    end
end
end%Ends function heatOut