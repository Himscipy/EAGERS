function GenOutput = checkStartupCosts(GenOutput,StorPower,Alt,StartCost,dt,include,n)
%This function finds the n shortest segments of time the generator is on or
%off, and compares the start-up cost to the alternative configuration that
%avoids the start or re-start
% GenOutput is the bet dispatch at each time without considering start-up costs
% Alt is all of the other feasible combinations tested
% Binary is the current best on/off configuration
% StartCost is the startup cost of each generator
% dt is the duration of each time segment
% Type specifies the category of generator (currently this function is only interested in electric and CHP generators.
% include is what type of generators to consider (currently this is always electric and CHP generators)
% n is # of segments to check
global Plant
nS = length(Alt.Disp);
nG = length(Plant.Generator);
InitialOutput = GenOutput;
SkipOn = [];
SkipOff = [];
Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
Binary = true(nS+1,nG);
for t = 1:1:nS+1
    Binary(t,Dispatchable) = GenOutput(t,Dispatchable)>0;
end
inc = false(1,nG);
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,include)
        inc(i) = true;
    end
end
I = ones(length(dt),1);
for t = 1:1:nS
    B = Binary(t+1,inc);
    i = find(ismember(Alt.Binary{t}(:,inc),B,'rows'));
    if length(i)>1
        I(t) = i(1);
        disp('Multiple binary options in checkStartupCosts')
    elseif ~isempty(i)
        I(t) = i;
    else
        disp('No binary options in checkStartupCosts')
    end
end
%find I(t), index that results in lowest cost once start-up is considered.
seg = 1;
onSeg = 1;
while ~isempty(onSeg) && seg<n
    %%Try and remove the shortest generator on segment
    [onSeg, offSeg] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff,inc);
    %% need to sort segments by length, StartCost and operating cost
    
    if isempty(onSeg) && isempty(offSeg)
        break %no  segments to check
    elseif ~isempty(onSeg)
        [~,Ion] = min(onSeg(:,4)-StartCost(onSeg(:,1))'/(2*max(StartCost))+onSeg(:,2)/nS); %shortest segment, and if equal lengths, highest start cost, if same generator then first segment
        k = onSeg(Ion,1); %generator index
        t1 = onSeg(Ion,2); %index when generator turns on
        t2 = onSeg(Ion,3);%index when generator shuts off
        %% Find the cheapest feasible alternative dispatch without this generator (only use generators that were on previously or will be on)
        %First try alternate generators at same time step that may have cheaper startup (or that are already on)
        [I,Binary,NoAlt] = altGeneration(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,dt);
        if NoAlt
            %Second, can it get to the final storage capacity without this on-segment?
            [Alt,I,Binary,mustReplace] = avoidGeneration(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,dt);
            if mustReplace
                %Third try moving generation earlier or later with same generator (if there is storage)
                [Alt,I,Binary,cantMove] = moveGeneration(I,k,t1,t2,Alt,GenOutput,Binary,inc,dt);
                if cantMove
                    SkipOn(end+1,1:4) = onSeg(Ion,:); %add to list of segments to avoid
                end
            end
        end
        GenOutput = updateStorageState(GenOutput,Alt,I,StorPower,dt);
    end
    %%Try and remove the shortest generator off segment, if it is a shorter segment than the next shortest on segment
    [onSeg, offSeg] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff,inc);
    if isempty(onSeg) && isempty(offSeg)
        break %no  segments to check
    elseif ~isempty(offSeg) && (isempty(onSeg) || (min(offSeg(:,4)-StartCost(offSeg(:,1))'/(2*max(StartCost))+offSeg(:,2)/nS) < min(onSeg(:,4)-StartCost(onSeg(:,1))'/(2*max(StartCost))+onSeg(:,2)/nS)))%preference is to turn things off, rather than turn things on
        %third leave generator on, maybe turn off a different generator
        [~,Ioff] = min(offSeg(:,4)-StartCost(offSeg(:,1))'/(2*max(StartCost))+offSeg(:,2)/nS); %shortest segment, and if equal lengths, highest start cost
        k = offSeg(Ioff,1); %generator index
        t1 = offSeg(Ioff,2); %index when generator turns on
        t2 = offSeg(Ioff,3);%index when generator shuts off
        [I,Binary,CantKeepOn] = leaveGenOn(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,inc,dt);
        if CantKeepOn
            SkipOff(end+1,1:4) = offSeg(Ioff,:); %add to list of segments to avoid
        end
        GenOutput = updateStorageState(GenOutput,Alt,I,StorPower,dt);
    end 
    seg = seg+1;
end
diff = GenOutput - InitialOutput;
end% Ends function checkStartupCosts

function GenOutput = updateStorageState(GenOutput,Alt,I,StorPower,dt)
%pull the corresponding best dispatches with the start-cost considered
global Plant
status = GenOutput(1,:);
nG = length(Plant.Generator);
nS = length(Alt.Disp);
for t = 1:1:nS
    newStatus = Alt.Disp{t}(I(t),:);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
            loss = (Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize*Plant.Generator(i).QPform.Stor.DischEff);
            if StorPower(t,i)>0 %discharging
                d_SOC = (-StorPower(t,i)/Plant.Generator(i).QPform.Stor.DischEff - loss)*dt(t);
            else %charging
                d_SOC = (-StorPower(t,i)*Plant.Generator(i).QPform.Stor.ChargeEff - loss)*dt(t);
            end
            new_d_SOC = -newStatus(i)*dt(t)/Plant.Generator(i).QPform.Stor.DischEff;
            newStatus(i) = status(i) + d_SOC + new_d_SOC;
            if newStatus(i)>Plant.Generator(i).QPform.Stor.UsableSize
                %Changed binary, but there was slack in other timesteps so storage will not actually overcharge
%                     disp(strcat('Warning, Potentially Over Charging Storage in ',num2str(sum(dt(1:t))),'_hours'))
                newStatus(i) = Plant.Generator(i).QPform.Stor.UsableSize;
            elseif newStatus(i)<0
                %Changed binary, but there was spare capacity in other timesteps so storage will not actually deplete
%                     disp(strcat('Warning, Potentially Depleating Storage in ',num2str(sum(dt(1:t))),'_hours'))
                newStatus(i) = 0;
            end
        end
    end
    status = newStatus;
    GenOutput(t+1,:) = status;
end
end%Ends function updateStorageState

function [Alt,I,Binary,mustReplace] = avoidGeneration(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,dt)
%need to re-do this for variable time steps
%looking for enough slack in other generators to avoid the overdepleating storage at the minimum, and to get to the same final state
global Plant DateSim
nS = length(Alt.Disp);
nG = length(GenOutput(1,:));
mustReplace = true;
inc = false(1,nG);
if strcmp(Plant.Generator(k).Type,'Chiller')
    include = {'Chiller'};
    out = 'C';
    output = 'Cooling';
elseif strcmp(Plant.Generator(k).Type,'Heater')
    include = {'Heater'};
    out = 'H';
    output = 'Heating';
elseif ismember(Plant.Generator(k).Type,{'CHP Generator';'Electric Generator';})
    include = {'CHP Generator';'Electric Generator';};
    out = 'E';
    output = 'Electricity';
end
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,include) 
        inc(i) = true;
    end
end
[UsefulStoredEnergy,~,stor] = StorState(GenOutput,dt);
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
Date = [DateSim; DateSim+Time/24];
dt = Time - [0; Time(1:end-1)];
scaleCost = updateGeneratorCost(Date); 
[~,~,spareGenCumulative] = GenLimit(GenOutput,Binary,dt);
rmvCost = 0;
if any(UsefulStoredEnergy.(out))>0
    MarginCost = MarginalCapacityCost(GenOutput,Date);
    remGen = zeros(nS,1);
    remHeat = zeros(nS,1);
    for t = t1:1:t2-1
        remGen(t:end) = remGen(t:end) + GenOutput(t+1,k)*dt(t);%need to replace this energy in the storage by the end, and replace enough early so that UsefulStoredEnergy - remStor + makeup does not go negative
        rmvCost = rmvCost + scaleCost(t+1,k)*GenOutput(t+1,k)./interp1(Plant.Generator(k).Output.Capacity*Plant.Generator(k).Size,Plant.Generator(k).Output.(output),GenOutput(t+1,k))*dt(t);
        if isfield(Plant.Generator(k).QPform,'constCost')
            rmvCost = rmvCost + Plant.Generator(k).QPform.constCost*scaleCost(t+1,k)*dt(t);
        end
        if strcmp(Plant.Generator(k).Type,'Chiller') && isfield(Plant.Generator(k).QPform.constDemand,'E')
            rmvCost = rmvCost + Plant.Generator(k).QPform.constDemand.E*min(nonzeros(MarginCost.E.Cost.SpinReserve(:,t,1)))*dt(t);
        end
        if strcmp(Plant.Generator(k).Type,'CHP Generator')
            remHeat(t:end) = remHeat(t:end) + GenOutput(t+1,k)*dt(t);
            D = 0;
            j = 1;
            if  GenOutput(t+1,k)>0
                remHeat(t:end) = remHeat(t:end) - Plant.Generator(k).QPform.constDemand.H;      
            end
            states = Plant.Generator(k).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(k).QPform.states(:,end))),end);
            while j<=length(states)
                d = min(GenOutput(t+1,k) - D,Plant.Generator(k).QPform.(states{j}).ub(2));
                D = D + d;
                remHeat(t:end) = remHeat(t:end) + d*Plant.Generator(k).QPform.output.H(j,2);
                j = j+1;
            end
        end
    end
    if all(remGen-spareGenCumulative.(out)<UsefulStoredEnergy.(out)) && spareGenCumulative.(out)(end)>=remGen(end) && (~strcmp(Plant.Generator(k).Type,'CHP Generator') || (all(remHeat-spareGenCumulative.H<UsefulStoredEnergy.H) && spareGenCumulative.H(end)>=remHeat(end)))
        if any(remGen>UsefulStoredEnergy.(out))%Need to make up generation before storage would go negative
            t_limit = min(nonzeros((1:nS)'.*(remGen>UsefulStoredEnergy.(out))));
        else
            t_limit = nS;
        end
        makeUpGen = 0;
        remE = max(remGen);
        sortMC = sortMarginalCost(MarginCost,out,t_limit,k,dt);
        if ~isempty(sortMC) && sortMC(end,1)>=remE
            makeUpCost = interp1(sortMC(:,1),sortMC(:,2),remE);
            newStorPower = zeros(nS,1);
            if (makeUpCost-rmvCost)<StartCost(k)%if spareGen is less than StartCost
                for t = t1:1:t2-1
                    Binary(t+1,k) = false; %turn off
                    newStorPower(t) = Alt.Disp{t}(I(t),k);
                    Alt.Disp{t}(I(t),k) = 0;%set output of this generator to zero
                end
                r = 1;
                while makeUpGen<remE && r<length(sortMC(:,1))
                    t = sortMC(r,4);
                    addGen = min(sortMC(r,5),remE);
                    Alt.Disp{t}(I(t),sortMC(r,3)) = Alt.Disp{t}(I(t),sortMC(r,3)) + addGen;
                    %%edit storage dispatch at this timestep so that UpdateStorageState works correctly
                    Alt.Disp{t}(I(t),:) = storAdd(Alt.Disp{t}(I(t),:),addGen,k);
                    makeUpGen = makeUpGen + addGen;
                    if t>=t1 && t<t2
                        newStorPower(t) = newStorPower(t) - addGen;
                    end
                    r = r+1;
                end
                for t = t1:1:t2-1
                    Alt.Disp{t}(I(t),stor.(out)) = Alt.Disp{t}(I(t),stor.(out))+newStorPower(t)/length(stor.(out));
                end
                mustReplace = false;
            end
        end
    end
end
end%Ends function avoidGeneration

function [Alt,I,Binary,cantMove] = moveGeneration(I,k,t1,t2,Alt,GenOutput,Binary,inc,dt)
%need to re-do this for variable time steps
global Plant
if ismember(Plant.Generator(k).Type,{'CHP Generator';'Electric Generator';})
    out = 'E';
elseif strcmp(Plant.Generator(k).Type,'Chiller')
    out = 'C';
elseif strcmp(Plant.Generator(k).Type,'Heater')
    out = 'H';
end
nS = length(Alt.Disp);
cantMove = true;
Ialt = I;
[UsefulStoredEnergy,~,~,SpareStorCap,~] = StorState(GenOutput,dt);
if any(UsefulStoredEnergy.(out))>0
    addStor = zeros(nS,1);
    remStor = zeros(nS,1);
    stops = nonzeros((1:nS)'.*(~Binary(2:end,k) & Binary(1:nS,k)));
    n = t2-t1; %number of steps generator is on for
    if any(stops<t1)
        %it was on previously, try to move earlier if storage permits
        t_stop = stops(find(stops<t1,1,'last'));
        %check if there is a feasible option to leave this generator on at earlier time steps
        for j = 0:1:n-1
            Opt = Alt.Binary{t_stop+j}; %all feasible options tested at this time
            BinaryNow = Opt(Ialt(t_stop+j),:);
            BinaryNow(k) = true;
            newBinary = nonzeros((1:length(Opt(:,1)))'.*ismember(Opt,BinaryNow,'rows'));
            if~isempty(newBinary)
                Ialt(t_stop+j) = newBinary;
            end
            addStor(t_stop+j:t1+j-1) = addStor(t_stop+j:t1+j-1) + GenOutput(t1+j+1,k);%need to hold this shifted energy in the storage from the previous time the generator shut down, until it had come on before
        end
        if ~any(Ialt(t_stop:t_stop+n-1) == I(t_stop:t_stop+n-1)) && all(addStor<SpareStorCap.(out))%possible to have generator on earlier
            I = Ialt; %use the alternative index
            for j = 0:1:n-1
                Binary(t_stop+j+1,inc) = Alt.Binary{t_stop+j}(I(t_stop+j),inc); %best alternative option tested at this time
                Alt.Disp{t_stop+j}(I(t_stop+j),k) = GenOutput(t1+j+1,k);
            end
            for t = t_stop+n:t2
                Binary(t+1,k) = false; %turn off
                Alt.Disp{t}(I(t),k) = 0;%set output of this generator to zero
            end
            cantMove = false;
        end
    end
    %if it cant move earlier, because of lack of storage space, try moving later
    if cantMove
        starts = nonzeros((1:nS)'.*(Binary(2:end,k) & ~Binary(1:nS,k)));
        t_start = starts(find(starts>t2,1,'first'))-n;%new starting point, n steps before its scheduled start
        if isempty(t_start)
            t_start = nS-n+1;%push as late as possible
        end
        %check if there is a feasible option to leave this generator on at earlier time steps
        for j = 0:1:n-1
            Opt = Alt.Binary{t_start+j}; %all feasible options tested at this time
            BinaryNow = Opt(Ialt(t_start+j),:);
            BinaryNow(k) = true;
            newBinary = nonzeros((1:length(Opt(:,1)))'.*ismember(Opt,BinaryNow,'rows'));
            if~isempty(newBinary)
                Ialt(t_start+j) = newBinary;
            end
            remStor(t1+j:t_start+j-1) = remStor(t1+j:t_start+j-1) + GenOutput(t1+j+1,k);%need to hold this shifted energy in the storage from the previous time the generator shut down, until it had come on before
        end
        if ~any(Ialt(t_start:t_start+n-1) == I(t_start:t_start+n-1)) && all(remStor>UsefulStoredEnergy.(out))%possible to have generator on later
            I = Ialt; %use the alternative index
            for j = 0:1:n-1
                Binary(t_start+j+1,inc) = Alt.Binary{t_start+j}(I(t_start+j),inc); %best alternative option tested at this time
                Alt.Disp{t_start+j}(I(t_start+j),k) = GenOutput(t1+j+1,k);
            end
            for t = t1:t_start-1
                Binary(t+1,k) = false; %turn 
                Alt.Disp{t}(I(t),k) = 0;%set output of this generator to zero
            end
            cantMove = false;
        end
    end
end
end%Ends function moveGeneration

function [I,Binary,NoAlt] = altGeneration(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,dt)
global Plant
NoAlt = true;
nS = length(Alt.Disp);
nG = length(StartCost);
inc = false(1,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(k).Type,'Chiller')
        out = 'C';
        include = {'Chiller'};
        if strcmp(Plant.Generator(i).Type,'Chiller') || (strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Cooling'))
            inc(i) = true;
        end
    end
    if strcmp(Plant.Generator(k).Type,'Heater')
        out = 'H';
        include = {'Heater'};
        if strcmp(Plant.Generator(i).Type,'Heater') || strcmp(Plant.Generator(i).Type,'CHP Generator')  || (strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat'))
            inc(i) = true;
        end
    end
    if strcmp(Plant.Generator(k).Type,'CHP Generator')
        out = 'E';
        include = {'CHP Generator';'Electric Generator';};
        dHeat = zeros(t2-t1,1);
        spareHeat = zeros(t2-t1,1);
        if strcmp(Plant.Generator(i).Type,'Heater') || strcmp(Plant.Generator(i).Type,'CHP Generator') || strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'Electric Storage') || (strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat'))
            inc(i) = true;
        end
    end
    if strcmp(Plant.Generator(k).Type,'Electric Generator')
        out = 'E';
        include = {'CHP Generator';'Electric Generator';};
        if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')  || strcmp(Plant.Generator(i).Type,'Electric Storage')
            inc(i) = true;
        end
    end
end
Ialt = I;
dCost = zeros(t2-t1,1);
dGen = zeros(t2-t1,1);
[UsefulStoredEnergy,StorGenAvail,~] = StorState(GenOutput,dt);
AltBinary = Binary;
for t = t1:1:t2-1
    Opt = Alt.Binary{t}; %all feasible options tested at this time
    cost2 = Alt.Cost{t}; %cost for these options
    cost2(Opt(:,k)) = inf;%make the cost infinite for all combinations with generator i
    cost2(ismember(Opt(:,inc),Binary(t+1,inc),'rows')) = inf;%make the cost infinite if its not changing the status of any of the included generator type
    for j = 1:1:nG
        if ~any(Binary(t1:t2+1,j))%would be adding a startup
            cost2 = cost2 + Opt(:,j)*StartCost(j);
        end
    end
    if any(~isinf(cost2))
        [dCost(t-t1+1),Ialt(t)] = min(cost2);
        dGen(t-t1+1) = Alt.Disp{t}(I(t),k) - Alt.Disp{t}(Ialt(t),k);
        AltBinary(t+1,:) = Alt.Binary{t}(Ialt(t),:);
    else
        dCost = inf;%no feasible combinations without turning on a different generator (don't make changes to this segment)
        break
    end
end
spareGen = zeros(t2-t1,1);
MaxOut = GenLimit(GenOutput,AltBinary,dt);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'CHP Generator')
        %convert GenOutput to heatOut
        heatOut = zeros(t2-t1,1);
        for t = t1:1:t2-1
            D = 0;
            j = 1;
            if  GenOutput(t+1,i)>0
                heatOut(t-t1+1) = heatOut(t-t1+1) - Plant.Generator(i).QPform.constDemand.H;      
            end
            states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
            while j<=length(states)
                d = min(GenOutput(t+1,i) - D,Plant.Generator(i).QPform.(states{j}).ub(2));
                D = D + d;
                heatOut(t-t1+1) = heatOut(t-t1+1) + d*Plant.Generator(i).QPform.output.H(j,2);
                j = j+1;
            end
        end
    end
    if i~=k && ismember(Plant.Generator(i).Type,include) && Alt.Disp{t}(Ialt(t),i)>0
        for t = t1:1:t2-1
            spareGen(t-t1+1) = spareGen(t-t1+1) + MaxOut.(out)(t+1,i) - GenOutput(t+1,i);%capacity is either UB or limited by ramping
        end
    elseif strcmp(Plant.Generator(k).Type,'Heater') && strcmp(Plant.Generator(i).Type,'CHP Generator')
        for t = t1:1:t2-1
            spareGen(t-t1+1) = spareGen(t-t1+1) + MaxOut.H(t+1,i) - heatOut(t-t1+1);
        end
    end
    if strcmp(Plant.Generator(k).Type,'CHP Generator')
        if i==k
            dHeat(t-t1+1) = heatOut(t-t1+1) - 0;
        end
        if strcmp(Plant.Generator(i).Type,'Heater')
            for t = t1:1:t2-1
                spareHeat(t-t1+1) = spareHeat(t-t1+1) + MaxOut.H(t+1,i) - GenOutput(t+1,i);
            end
        elseif strcmp(Plant.Generator(i).Type,'CHP Generator')
            for t = t1:1:t2-1
                spareHeat(t-t1+1) = spareHeat(t-t1+1) + MaxOut.H(t+1,i) - heatOut(t-t1+1);
            end
        end
    end
end
if sum(dCost)<StartCost(k) && all(dGen<(spareGen+StorGenAvail.(out)(t1:t2-1)./dt(t1:t2-1))) && sum(dGen-spareGen)<min(UsefulStoredEnergy.(out)(t1:end)) %sum of the marginal increase in cost is less than the start-up cost, there is spare capacity in the other generators & storage, and the cumulative loss of generation does not deplete the storage
    if ~strcmp(Plant.Generator(k).Type,'CHP Generator') || (all(dHeat<(spareHeat+StorGenAvail.H(t1:t2-1)./dt(t1:t2-1))) && sum(dHeat-spareHeat)<min(UsefulStoredEnergy.H(t1:end)))
        I = Ialt; %use the alternative index
        for t = t1:1:t2-1
            Binary(t+1,inc) = Alt.Binary{t}(I(t),inc); %best alternative option tested at this time
        end
        NoAlt = false;
    end
end
end%Ends function altGeneration

function [I,Binary,CantKeepOn] = leaveGenOn(I,k,t1,t2,StartCost,Alt,GenOutput,Binary,inc,dt)
%% Find the cheapest feasible alternative dispatch that keeps this generator on (only use generators that were on previously or will be on)
global Plant
if ismember(Plant.Generator(k).Type,{'CHP Generator';'Electric Generator';})
    out = 'E';
elseif strcmp(Plant.Generator(k).Type,'Chiller')
    out = 'C';
elseif strcmp(Plant.Generator(k).Type,'Heater')
    out = 'H';
end
CantKeepOn = true;
nS = length(Alt.Disp);
nG = length(inc);
%only allow generators that are on at begining or end to be involved (or that have smaller start-up cost)
Ialt = I;
dCost = zeros(t2-t1,1);
dGen = zeros(t2-t1,1);
slackGen = zeros(t2-t1,1);
[~,~,~,SpareStorCap,StorSlackAvail] = StorState(GenOutput,dt);
if t1 == 1
    Prev = GenOutput(1,:);
else
    Prev = Alt.Disp{t1-1}(I(t1-1),:);
end
for t = t1:1:t2-1
    Opt = Alt.Binary{t}; %all feasible options tested at this time
    cost2 = Alt.Cost{t}; %cost for these options
    cost2(~Opt(:,k)) = inf;%make the cost infinite for all combinations without generator i
    cost2(ismember(Opt(:,inc),Binary(t+1,inc),'rows')) = inf;%make the cost infinite if its not changing the status of any of the included generator type
    for j = 1:1:nG
        if ~any(Binary(t1:t2+1,j))%would be adding a startup
            cost2 = cost2 + Opt(:,j)*StartCost(j);
        end
    end
    if any(~isinf(cost2))
        [dCost(t-t1+1),Ialt(t)] = min(cost2);
        [dGen(t-t1+1),slackGen(t-t1+1)] = slackCap(Alt.Disp{t}(I(t),:),Alt.Disp{t}(Ialt(t),:),Prev,k,sum(dt(t1:t)));
    else
        dCost = inf;%no feasible combinations keeping this generator active (don't make changes to this segment)
        break
    end
end
if sum(dCost)<StartCost(k) && all(dGen<(slackGen+StorSlackAvail.(out)(t1:t2-1)./dt(t1:t2-1))) && sum(dGen-slackGen)<min(SpareStorCap.(out)(t1:end)) %sum of the marginal increase in cost is less than the start-up cost, there is spare capacity in the other generators & storage, and the cumulative loss of generation does not deplete the storage
    I = Ialt; %use the alternative index
    for t = t1:1:t2-1
        Binary(t+1,:) = Alt.Binary{t}(I(t),:); %best alternative option tested at this time
    end
    CantKeepOn = false;
end
end%Ends function leaveGenOn

function [onSeg, offSeg] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff,inc)
onSeg = [];
offSeg = [];
nG = length(StartCost);
nS = length(dt);
Horizon = sum(dt);
% find length of segments that a generator is on or off
for i = 1:1:nG
    if StartCost(i)>0 && any(~Binary(:,i)) && inc(i)
        starts = nonzeros((1:nS)'.*(Binary(2:end,i) & ~Binary(1:nS,i)));
        stops = nonzeros((1:nS)'.*(~Binary(2:end,i) & Binary(1:nS,i)));
        nOn = length(starts);
        nOff = length(stops);
        if nOn>0 && nOff>0 %only look at generators that turn both on and off during the window
            if stops(1)<starts(1) %generator is off for a segment
                Seg = [i, stops(1), starts(1), sum(dt(stops(1):starts(1)-1))];
                if (isempty(SkipOff) || ~any(ismember(SkipOff,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    offSeg(end+1,1:4) = Seg;
                end
                stops = stops(2:end);
                nOff  = nOff - 1;
            end
            j = 0;
            while j < nOff
                j = j + 1;
                Seg = [i, starts(j), stops(j), sum(dt(starts(j):stops(j)-1))]; %index of generator, start index, stop index, duration of segment
                if (isempty(SkipOn) || ~any(ismember(SkipOn,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    onSeg(end+1,1:4) = Seg;
                end
                if j<nOn
                    Seg = [i, stops(j), starts(j+1), sum(dt(stops(j):starts(j+1)-1))];
                    if (isempty(SkipOff) || ~any(ismember(SkipOff,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                        offSeg(end+1,1:4) = Seg;
                    end
                end
            end
        end
    end
end
end%Ends function segmentLength


function [dGen,slackGen] = slackCap(Original,New,Prev,k,dt)
global Plant
nG = length(Plant.Generator);
dGen = Plant.Generator(k).QPform.A.lb(end);
if strcmp(Plant.Generator(k).Type,'Chiller')
    include = {'Chiller'};
elseif strcmp(Plant.Generator(k).Type,'Heater')
    include = {'Heater'};
elseif ismember(Plant.Generator(k).Type,{'CHP Generator';'Electric Generator';})
    include = {'CHP Generator';'Electric Generator';};
end
slackGen = 0;
for i = 1:1:nG
    if i~=k && ismember(Plant.Generator(i).Type,include) 
        if Original(i)>0 && New(i) >0
            slackGen = slackGen + min(Original(i)-Plant.Generator(i).QPform.A.lb(end),Original(i)-Prev(i)+Plant.Generator(i).VariableStruct.dX_dt*dt);%slack capacity is either Original - LB or limited by ramping
        elseif Original(i)>0
            slackGen = slackGen + Original(i);
        elseif New(i)>0
            slackGen = slackGen - New(i);
        end
    end
end
end%Ends function slackCap

function [StoredEnergy,StorGenAvail,stor,SpareStorCap,StorSlackAvail] = StorState(GenOutput,dt)
global Plant
nG = length(Plant.Generator);
nS = length(GenOutput(2:end,1));
StoredEnergy.E = zeros(nS,1);
StoredEnergy.H = zeros(nS,1);
StoredEnergy.C = zeros(nS,1);
StorGenAvail = StoredEnergy;
SpareStorCap = StoredEnergy;
StorSlackAvail = StoredEnergy;
stor.E = [];
stor.H = [];
stor.C = [];
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Cooling')
        buff = Plant.Generator(i).QPform.Stor.UsableSize*(Plant.Generator(i).VariableStruct.Buffer/100);
        StoredEnergy.C = StoredEnergy.C + (GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff;
        StorGenAvail.C = StorGenAvail.C +  min((GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff./dt,Plant.Generator(i).QPform.Stor.PeakDisch);
        SpareStorCap.C = SpareStorCap.C + (Plant.Generator(i).QPform.Stor.UsableSize - buff - GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff;
        StorSlackAvail.C = StorSlackAvail.C +  min((Plant.Generator(i).QPform.Stor.UsableSize-buff-GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff./dt,Plant.Generator(i).QPform.Stor.PeakCharge);
        stor.C(end+1) = i;
    end
    if strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat')
        buff = Plant.Generator(i).QPform.Stor.UsableSize*(Plant.Generator(i).VariableStruct.Buffer/100);
        StoredEnergy.H = StoredEnergy.H + (GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff;
        StorGenAvail.H = StorGenAvail.H +  min((GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff./dt,Plant.Generator(i).QPform.Stor.PeakDisch);
        SpareStorCap.H = SpareStorCap.H + (Plant.Generator(i).QPform.Stor.UsableSize - buff - GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff;
        StorSlackAvail.H = StorSlackAvail.H +  min((Plant.Generator(i).QPform.Stor.UsableSize-buff-GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff./dt,Plant.Generator(i).QPform.Stor.PeakCharge);
        stor.H(end+1) = i;
    end
    if  strcmp(Plant.Generator(i).Type,'Electric Storage')
        buff = Plant.Generator(i).QPform.Stor.UsableSize*(Plant.Generator(i).VariableStruct.Buffer/100);
        StoredEnergy.E = StoredEnergy.E + (GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff;
        StorGenAvail.E = StorGenAvail.E +  min((GenOutput(2:end,i)-buff)*Plant.Generator(i).QPform.Stor.DischEff./dt,Plant.Generator(i).QPform.Stor.PeakDisch);
        SpareStorCap.E = SpareStorCap.E + (Plant.Generator(i).QPform.Stor.UsableSize - buff - GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff;
        StorSlackAvail.E = StorSlackAvail.E +  min((Plant.Generator(i).QPform.Stor.UsableSize-buff-GenOutput(2:end,i))/Plant.Generator(i).QPform.Stor.ChargeEff./dt,Plant.Generator(i).QPform.Stor.PeakCharge);
        stor.E(end+1) = i;
    end
end
end%Ends function StorState


function Dispatch = storAdd(Prev,addGen,k)
global Plant
nG = length(Plant.Generator);
Dispatch = Prev;
for i = 1:1:nG
    if strcmp(Plant.Generator(k).Type,'Chiller') && strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Cooling')
        change = min(addGen, Plant.Generator(i).QPform.Stor.PeakCharge + Prev(i));
        addGen = addGen - change;
        Dispatch(i) = Prev(i) - change;
    end
    if strcmp(Plant.Generator(k).Type,'Heater') && strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat')
        change = min(addGen, Plant.Generator(i).QPform.Stor.PeakCharge + Prev(i));
        addGen = addGen - change;
        Dispatch(i) = Prev(i) - change;
    end
    if ismember(Plant.Generator(k).Type,{'CHP Generator';'Electric Generator';}) && strcmp(Plant.Generator(i).Type,'Electric Storage')
        change = min(addGen, Plant.Generator(i).QPform.Stor.PeakCharge + Prev(i));
        addGen = addGen - change;
        Dispatch(i) = Prev(i) - change;
    end
end
end%Ends function storAdd