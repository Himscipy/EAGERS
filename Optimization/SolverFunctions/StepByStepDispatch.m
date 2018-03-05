function GenOutput = StepByStepDispatch(Forecast,scaleCost,dt,FirstProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
% StorPower is the amount of power comming from (positive) or going into (negative) the storage device at each timestep according to the first dispatch
global Plant CurrentState
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
    T_build = CurrentState.Buildings;
else
    nB = 0;
    T_build = [];
end
QP.Organize.Enabled = ones(1,nG);%message of which components are not enabled
netDemand = agregateDemand(Forecast,nB);
if isempty(FirstProfile) %finding initial conditions
    limit = 'unconstrained';
    marginal = updateMarginalCost([],scaleCost,[],0);%update marginal cost
    QP = updateMatrices1Step(Plant.OneStep,Forecast,scaleCost,marginal,[],dt,[],[],limit,1,T_build);
    K = createCombinations(QP,netDemand,[],dt(1),1,[],'E');%% create a matrix of all possible combinations 
    [BestDispatch,~] = eliminateCombinations(QP,K,[],[],netDemand,dt(1),1);%% combination elimination loop
    [GenOutput,~] = updateIC([],BestDispatch,[],dt(1),limit);
    if isfield(Forecast,'Renewable')
        GenOutput(1,1:nG) = GenOutput(1,1:nG) + Forecast.Renewable;
    end
else
    limit = 'initially constrained';
    StartCost = zeros(1,nG);
    for i = 1:1:nG
        if ~Plant.Generator(i).Enabled
            QP.Organize.Enabled(i) = 0;
        end
        if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
            StartCost(i) = Plant.Generator(i).VariableStruct.StartCost;
        end
    end
    IC = CurrentState.Generators;%, CurrentState.Lines,CurrentState.Buildings(1,:)];
    [nS,~] = size(scaleCost);
    GenOutput = zeros(nS+1,length(IC));
    GenOutput(1,:) = IC;
    Alt.Disp = cell(nS,1);
    Alt.Cost = cell(nS,1);
    Alt.Binary = cell(nS,1);
    
    %If there are more than 2 chillers, solve it as a seperate problem
    [K_chill,Alt_C,FirstProfile] = seperateChillerProblem(Forecast,FirstProfile,netDemand,scaleCost,StartCost,IC,dt,T_build);
    StorPower = zeros(nS,nG);
    for t = 1:1:nS %for every timestep
        if isfield(netDemand,'H')
            vH = valueHeat([IC;FirstProfile(t+1,:)],netDemand.H(t),dt(t));
        else
            vH = valueHeat([IC;FirstProfile(t+1,:)],0,dt(t));
        end
        StorPower(t,:) = findStorPower([IC;FirstProfile(t+1,:)],dt(t));
        marginal = updateMarginalCost(FirstProfile(t+1,:),scaleCost(t,:),dt(t),vH);%update marginal cost
        QP = updateMatrices1Step(Plant.OneStep,Forecast,scaleCost(t,:),marginal,StorPower(t,:),dt,IC,FirstProfile,limit,t,T_build);
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
                out = fieldnames(Plant.Generator(i).QPform.output);
                if isfield(netDemand,out{1})
                    netDemand.(out{1})(t) = netDemand.(out{1})(t) - StorPower(t,i);
                end
            end
        end
        EC = FirstProfile(t+1,:);
        Prev = GenOutput(t,:);
        K = createCombinations(QP,netDemand,EC,dt(t),t,K_chill(t+1,:),'E');%% create a matrix of all possible combinations using the best chiller dispatch combination
        [BestDispatch,Alt] = eliminateCombinations(QP,K,Alt,Prev,netDemand,dt(t),t);%% combination elimination loop
        
        %If there was not actually a feasible combination of generators with the best chiller combination, try the next best chiller combination
        j = 1;
        while isempty(Alt.Disp{t}) && ~isempty(Alt_C) && j<=length(Alt_C.Disp{t}(:,1))
            K = createCombinations(QP,netDemand,EC,dt(t),t,Alt_C.Binary{t}(j,:),'E');%% create a matrix of all possible combinations using the next best chiller dispatch combination
            [BestDispatch,Alt] = eliminateCombinations(QP,K,Alt,Prev,netDemand,dt(t),t);%% combination elimination loop
            j = j+1;
        end
        
        %update Initial conditions and building temperatures
        [GenOutput(t+1,:),IC] = updateIC(IC,BestDispatch,FirstProfile(t+1,:),dt(t),limit);%only updates the storage states in IC
        if ~isempty(T_build)
            Tzone = zeros(2,nB);
            Twall = zeros(2,nB);
            for i = 1:1:nB
                [Tzone(:,i),Twall(:,i)] = BuildingSimulate(Plant.Building(i),Forecast.Weather.Tdb(t),Forecast.Weather.RH(t),dt(t)*3600,Forecast.Building.InternalGains(t,i),Forecast.Building.ExternalGains(t,i),Alt.Cooling{t}(i),Alt.Heating{t}(i),Forecast.Building.AirFlow(t,i),Forecast.Building.Damper(t,i),T_build(1,i),T_build(2,i));
            end
            T_build = [Tzone(2,:);Twall(2,:)];
        end
    end
    
    %change best dispatch based on re-start costs
    include = {'Electric Generator';'CHP Generator';'Heater';'Electrolyzer';'Hydrogen Generator';'Cooling Tower';};
    if isempty(Alt_C)
        include(end+1) = {'Chiller'};
    end
    GenOutput = checkStartupCosts(GenOutput,StorPower,Alt,StartCost,dt,include,40);
    if isfield(Forecast,'Renewable')
        GenOutput(2:nS+1,1:nG) = GenOutput(2:nS+1,1:nG) + Forecast.Renewable;
    end
end
end% Ends function StepByStepDispatch

function [GenOutput,IC] = updateIC(IC,BestDispatch,FirstProfile,dt,limit)
global Plant
nG = length(Plant.Generator);
GenOutput = BestDispatch;
if strcmp(limit,'unconstrained')
    IC = [];
else
    EC = BestDispatch(1:nG);
    if ~isempty(FirstProfile)
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
                d_SOC = FirstProfile(i) - IC(i);%first profile already accounted for loss
                new_d_SOC = -BestDispatch(i)*dt/Plant.Generator(i).QPform.Stor.DischEff;
                EC(i) = IC(i) + d_SOC + new_d_SOC;
                if EC(i)<0
                    error = EC(i)/Plant.Generator(i).QPform.Stor.UsableSize*100;
                    EC(i) = 0;
                    disp(strcat('Warning: ',Plant.Generator(i).Name,'_ is going negative by_',num2str(error),'%'))
                elseif EC(i)>Plant.Generator(i).QPform.Stor.UsableSize
                    error = (EC(i)-Plant.Generator(i).QPform.Stor.UsableSize)/Plant.Generator(i).QPform.Stor.UsableSize*100;
                    EC(i) = Plant.Generator(i).QPform.Stor.UsableSize;
                    disp(strcat('Warning: ',Plant.Generator(i).Name,'_ is exceeding max charge by_',num2str(error),'%'))
                end
            end
        end
    end
    GenOutput(1,1:nG) = EC;
    GenOutput(1,nG+1:end) = BestDispatch(nG+1:end);
    if strcmp(limit,'constrained')%if its constrained but not initially constrained then make the last output the initial condition
        IC = EC;
    else %if strcmp(limit,'initially constrained')
        for i = 1:1:nG
            if isfield(Plant.Generator(i).QPform,'Stor')
                IC(i) = EC(i); %update the state of storage
            end
        end
    end
end
end%ends function updateIC

function [K_chill,Alt_C,FirstProfile] = seperateChillerProblem(Forecast,FirstProfile,netDemand,scaleCost,StartCost,IC,dt,T_build)
%%optimize chilling dispatch first
global Plant
Alt_C = [];
nC = 0;
Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
nG = length(Plant.Generator);
for i = 1:1:nG
    if Plant.Generator(i).Enabled && Dispatchable(i) 
        if ismember(Plant.Generator(i).Type,{'Chiller'})
            nC = nC+1;
        end
    end
end
[nS,~] = size(scaleCost);
K_chill = zeros(nS+1,0);
StorPower = zeros(nS,nG);
if nC>2 && any(netDemand.C>0)%more than 2 dispatchable chillers, perform seperate optimization and keep n_test best scenarios
    GenOutput = zeros(nS+1,length(IC));
    GenOutput(1,:) = IC;
    K_chill = zeros(nS+1,nG);
    Alt_C.Disp = cell(nS,1);
    Alt_C.Cost = cell(nS,1);
    Alt_C.Binary = cell(nS,1);
    D.C = netDemand.C;
    for t = 1:1:nS %for every timestep
        if isfield(netDemand,'H')
            vH = valueHeat([IC;FirstProfile(t+1,:)],netDemand.H(t),dt(t));
        else
            vH = valueHeat([IC;FirstProfile(t+1,:)],0,dt(t));
        end
        StorPower(t,:) = findStorPower([IC;FirstProfile(t+1,:)],dt(t));
        marginal = updateMarginalCost(FirstProfile(t+1,:),scaleCost(t,:),dt(t),vH);%update marginal cost
        QP = updateMatrices1Step(Plant.OneStep,Forecast,scaleCost(t,:),marginal,StorPower(t,:),dt,IC,FirstProfile,'initially constrained',t,T_build);
        QP_C = chillerOnly1step(QP,marginal,FirstProfile(t+1,:));
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Thermal Storage';}) && isfield(Plant.Generator(i).QPform.output,'C')
                D.C(t) = D.C(t) - StorPower(t,i);
            end
        end
        K_C = createCombinations(QP_C,D,FirstProfile(t+1,:),dt(t),t,[],'C');%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
        [BestDispatch,Alt_C] = eliminateCombinations(QP_C,K_C,Alt_C,IC,D,dt(t),t);%% combination elimination loop
        [GenOutput(t+1,:),IC] = updateIC(IC,BestDispatch,FirstProfile(t+1,:),dt(t),'initially constrained');%only updates the storage states in IC
    end
    for i = 1:1:nG
        if ~(strcmp(Plant.Generator(i).Type,'Chiller') || ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';}))
            GenOutput(2:end,i) = FirstProfile(2:end,i);
        end
    end
    GenOutput = checkStartupCosts(GenOutput,StorPower,Alt_C,StartCost,dt,{'Chiller';},40);
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'Chiller') ||((ismember(Plant.Generator(i).Type,{'Thermal Storage';}) && isfield(Plant.Generator(i).QPform.output,'C')))
            FirstProfile(2:end,i) = GenOutput(2:end,i);
        end
    end
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'Thermal Storage';'Utility';'Electric Storage';'Solar';'Wind';})
            K_chill(:,i) = 1;  
        else
            K_chill(:,i) = FirstProfile(:,i)>0;
        end
    end
end
end %Ends function seperateChillerProblem

function netDemand = agregateDemand(Forecast,nB)
netDemand = [];
if isfield(Forecast,'Demand')
    Outs = fieldnames(Forecast.Demand);
    for j = 1:1:length(Outs)
        netDemand.(Outs{j}) = sum(Forecast.Demand.(Outs{j}),2);
    end
end
for i = 1:1:nB
    Outs2 = {'E';'H';'C';};
    nS = length(Forecast.Timestamp);
    for j = 1:1:length(Outs2)
        if ~isfield(netDemand,Outs2{j})
            netDemand.(Outs2{j}) = zeros(nS,1);
        end
        netDemand.(Outs2{j}) = netDemand.(Outs2{j}) + Forecast.Building.(strcat(Outs2{j},'0'))(:,i);
    end
end
end%Ends function agregateDemand

function vH = valueHeat(Dispatch,excessHeat,dt)
global Plant
nG = length(Plant.Generator);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat')
        loss = dt*(Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
        d_SOC = Dispatch(2,i) - Dispatch(1,i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
        if (d_SOC/dt + loss)<0 %discharging
            excessHeat = excessHeat -(d_SOC/dt + loss)*Plant.Generator(i).QPform.Stor.DischEff;
        else %charging
            excessHeat = excessHeat -(d_SOC/dt +loss)/Plant.Generator(i).QPform.Stor.ChargeEff; 
        end
    end
    if strcmp(Plant.Generator(i).Type,'CHP Generator')
        j = 1;
        D = 0;
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        if Dispatch(2,i)>0
            excessHeat = excessHeat - Plant.Generator(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
        end
        while j<=length(states) && Dispatch(2,i)-D>0
            x_i = min(Dispatch(2,i)-D,Plant.Generator(i).QPform.(states{j}).ub(end));
            excessHeat = excessHeat + x_i*Plant.Generator(i).QPform.output.H(min(j,length(Plant.Generator(i).QPform.output.H(:,1))),1);
            D = D + x_i;
            j = j+1;
        end
    end
end
vH = excessHeat<=1e-1;
end%Ends function valueHeat

function StorPower = findStorPower(Dispatch,dt)
global Plant
nG = length(Plant.Generator);
StorPower = zeros(length(dt),nG);  
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Stor')
        for t = 1:1:length(dt)
            loss = dt(t)*(Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
            d_SOC = Dispatch(t+1,i) - Dispatch(t,i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
            if (d_SOC/dt(t) + loss)<0 %discharging
                StorPower(t,i) = -(d_SOC/dt(t) + loss)*Plant.Generator(i).QPform.Stor.DischEff;
            else %charging
                StorPower(t,i) = -(d_SOC/dt(t) +loss)/Plant.Generator(i).QPform.Stor.ChargeEff; 
            end
        end
    end
end
end%Ends function findStorPower