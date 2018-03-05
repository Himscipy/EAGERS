function K = createCombinations(QP,netDemand,IC,dt,t,K_C,Out)
%% test if each combination is capable of meeting demand
%if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
%if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
%%could improve this section to check if feasible with network losses (not sure how)
global Plant 
if strcmp(Out,'C')
    C = netDemand.C;
    netDemand = [];
    netDemand.C = C;
    include = {'Chiller'};
    req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    [K,~] = buildK(include,QP.Organize.Dispatchable,[]);
    [K,~] = UpperLowerLimit(K,QP,req,netDemand.C(t),[],IC,dt);
elseif strcmp(Out,'E')
    if isfield(netDemand,'E') 
        netDemand.E = netDemand.E;
    else 
        netDemand.E = 0;
    end
    if isfield(netDemand,'DC')
        netDemand.E = netDemand.E + netDemand.DC;
    end
    include = {'Electric Generator';'CHP Generator';'Hydrogen Generator';'Heater';'Electrolyzer'};%no utilities or storage systems are dispatchable
    if isempty(K_C)
        include(end+1:end+3) = {'Chiller';'Absorption Chiller';'Cooling Tower';};
    else
        netDemand = rmfield(netDemand,'C');
    end
    req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    
    nG = length(Plant.Generator);
    isCHP = false;
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'CHP Generator')
            isCHP = true;
        end
    end
    [K_all,~] = buildK(include,QP.Organize.Dispatchable,[]);
    [K_all,~] = UpperLowerLimit(K_all,QP,req,netDemand.E(t),[],IC,dt);
    if isempty(K_C)
       K = K_all;
       if isCHP && ~isempty(K)
            req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            [K,~] = UpperLowerLimit(K,QP,req,netDemand.H(t),'H',IC,dt);
       end
       if isfield(netDemand,'C') && ~isempty(K)
           req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cooling demand
           [K,~] = UpperLowerLimit(K,QP,req,netDemand.C(t),'C',IC,dt);
       end
    else
        [K,inc] = buildK(include,QP.Organize.Dispatchable,K_C);
        [K,notFeas] = UpperLowerLimit(K,QP,req,netDemand.E(t),[],IC,dt);
        if isempty(K) && ~isempty(notFeas)
            for i = 1:1:length(notFeas(:,1))%%need extra logic in this create combinations so that a set of power/heat generators that is not feasible with the ideal chiller set may be feasible with a different chiller set that also meets the cooling constraints
                if any(notFeas(i,inc'))
                    f_rows = ismember(K_all(:,inc'),notFeas(i,inc'),'rows');%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                else
                    f_rows = all(~K_all(:,inc'),2);%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                end
                K(end+1:end+nnz(f_rows),:) = K_all(f_rows,:);%add these rows to K
            end
        end
        if isCHP && ~isempty(K)
            req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            [K,notFeas] = UpperLowerLimit(K,QP,req,netDemand.H(t),'H',IC,dt);
            if isempty(K) && ~isempty(notFeas) && ~isempty(K_all)
                for i = 1:1:length(notFeas(:,1))
                    if any(notFeas(i,inc'))
                        f_rows = ismember(K_all(:,inc'),notFeas(i,inc'),'rows');%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                    else
                        f_rows = all(~K_all(:,inc'),2);%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                    end
                    K(end+1:end+nnz(f_rows),:) = K_all(f_rows,:);%add these rows to K
                end
            end
        end
    end
end
if isempty(K)
    disp('Zero feasible combinations to test: ERROR')
end
%All combinations at this point should be feasible

% Premise (Maybe?), any extra generators beyond what is necessary to meet demand, or
% that are more expensive than grid are unneccessary
%if a combination is feasible, and another combination includes all those
%generators + additional  generators, those other cases can be removed
end%Ends function createCombinations

function [K,notFeas] = UpperLowerLimit(K,QP,req,netDemand,Out,IC,dt)
global Plant
nG = length(K(1,:));
limitL = zeros(1,nG);
limitU = zeros(1,length(limitL));
for i = 1:1:length(limitL)
    states = QP.Organize.States{i};
    if ~isempty(states) && any(any(QP.Aeq(req,states)~=0))
        I = find(any(QP.Aeq(req,states)~=0));
        if strcmp(Plant.Generator(i).Type,'Utility')%if it is a utility
            if length(states) ==2
                limitL(i) = sum(QP.Aeq(req(I),states(2))*QP.ub(states(2)));%if there is sellback, take 2nd state
                limitU(i) = sum(QP.Aeq(req(I),states(1))*QP.ub(states(1)));
            else
                limitL(i) = sum(QP.Aeq(req(I),states).*QP.lb(states)');
                limitU(i) = sum(QP.Aeq(req(I),states).*QP.ub(states)');
            end
        elseif any(strcmp(Plant.Generator(i).Type,{'Hydro Storage';'Electric Storage';'Thermal Storage';}))%if it is storage
            s = states(1);
            limitL(i) = sum(QP.Aeq(req(I),s)*QP.lb(s));
            limitU(i) = sum(QP.Aeq(req(I),s)*QP.ub(s));
            if ~isempty(IC)%bounds reduced by state of charge
                limitU(i) = min(IC(i).*QP.Aeq(req,s)/dt,limitU(i));
                chargingSpace = sum((Plant.Generator(i).QPform.Stor.UsableSize-IC(i)).*QP.Aeq(req,s));
                limitL(i) = max(-chargingSpace/dt, limitL(i));
            end
        else
            limitL(i) = sum(QP.Aeq(req(I),states).*QP.lb(states)');
            limitU(i) = sum(QP.Aeq(req(I),states).*QP.ub(states)');
            if limitU(i)<limitL(i)%if it consumes this demand (i.e. a chiller when examining electrical demand)
                temp = limitU(i);
                limitU(i) = limitL(i);
                limitL(i) = temp;
            end
        end
        if isfield(Plant.Generator(i).QPform,'constDemand') && ~isempty(Out) && isfield(Plant.Generator(i).QPform.constDemand,Out)
            limitU(i) = limitU(i) - Plant.Generator(i).QPform.constDemand.(Out);
            limitL(i) = limitL(i) - Plant.Generator(i).QPform.constDemand.(Out);
        end
    end
end
if ~isempty(Out) && strcmp(Out,'H') && Plant.optimoptions.excessHeat
    limitL = limitL - inf;
end
n = length(K(:,1));
A = (ones(n,1)*limitU);
A(K==0) = 0;
sumUB = sum(A,2);
B = (ones(n,1)*limitL);
B(K==0) = 0;
sumLB = sum(B,2);

% %% Trying to see if we can add transmission losses
% if ~isempty(QP.Organize.Transmission)
%     [sumUB,sumLB] = TransferLoss(QP,A,B,sumUB,sumLB);
% end 

keep = ((sumUB>=netDemand) & (sumLB<=netDemand));%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
notFeas = K(~keep,:);
K = K(keep,:);
end%ends function UpperLowerLimit

function [K,inc] = buildK(include,Dispatchable,K_C)
global Plant
nG = length(Plant.Generator);
inc = false(nG,1);
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,include) && Dispatchable(i)
        inc(i) = true;
    end
end
ninc = nonzeros((1:nG)'.*(~inc));
inc = nonzeros((1:nG)'.*inc);
n = length(inc);
K = zeros(2^n,nG); %K is a matrix of all the possible generator on/off combinations 

if ~isempty(ninc)
    for j = 1:1:length(ninc)
        if ~isempty(K_C) && strcmp(Plant.Generator(ninc(j)).Type,'Chiller')
            if K_C(ninc(j))
                K(:,ninc(j)) = ninc(j); % Chiller status determined earlier
            end
        else
            K(:,ninc(j)) = ninc(j); % all systems that are always included
        end
    end
end

for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=2^n
        K(r+1:r+z,inc(j)) = inc(j);
        r=r+2*z;
    end
end
end%ends function buildK