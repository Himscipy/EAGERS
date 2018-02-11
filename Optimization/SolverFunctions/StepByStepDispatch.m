function GenOutput = StepByStepDispatch(Forecast,scaleCost,dt,limit,FirstProfile)
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
nL = length(Plant.OneStep.Organize.States)-nG-nB;
Alt = [];
Alt_C = [];

Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
if ~isempty(FirstProfile)
    IC = [CurrentState.Generators, CurrentState.Lines,CurrentState.Buildings(1,:)];
    [nS,~] = size(scaleCost);
    StorPower = zeros(nS,nG);
    GenOutput = zeros(nS+1, nG+nL+nB);%nS should equal 1 in this case (finding IC)
    GenOutput(1,:) = IC;
    Alt.Disp = cell(nS,1);
    Alt.Cost = cell(nS,1);
    Alt.Binary = cell(nS,1);
else
    nS = 1;
    IC = [];
end
StartCost = zeros(1,nG);
QP.Organize.Enabled = zeros(1,nG);%message of which components are not enabled
nC = 0;
nH = 0;
for i = 1:1:nG
    if Plant.Generator(i).Enabled
        QP.Organize.Enabled(i) = 1;
        if Dispatchable(i) 
            if ismember(Plant.Generator(i).Type,{'Chiller'})
                nC = nC+1;
            elseif ismember(Plant.Generator(i).Type,{'Heater'})
                nH = nH+1;
            end
        end
    end
    if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
        StartCost(i) = Plant.Generator(i).VariableStruct.StartCost;
    end
end
netDemand = [];
if isfield(Forecast,'Demand')
    Outs = fieldnames(Forecast.Demand);
    for j = 1:1:length(Outs)
        netDemand.(Outs{j}) = sum(Forecast.Demand.(Outs{j}),2);
    end
end
for i = 1:1:nB
    Outs2 = {'E';'H';'C';};
    for j = 1:1:length(Outs2)
        if ~isfield(netDemand,Outs2{j})
            netDemand.(Outs2{j}) = zeros(nS,1);
        end
        netDemand.(Outs2{j}) = netDemand.(Outs2{j}) + Forecast.Building.(strcat(Outs2{j},'0'))(:,i);
    end
end
if isempty(FirstProfile) %finding initial conditions
    limit = 'unconstrained';
    marginal = instantMarginalCost([],scaleCost,[],[],[],1);%update marginal cost
    QP = updateMatrices1Step(Plant.OneStep,Forecast,scaleCost,marginal,[],dt,IC,[],limit,1,T_build);
end

%%optimize chilling dispatch first
K_chill = [];
if ~isempty(FirstProfile) && nC>2 && any(netDemand.C>0)%more than 2 dispatchable chillers, perform seperate optimization and keep n_test best scenarios
    K_chill = zeros(nS+1,nG);
    Alt_C.Disp = cell(nS,1);
    Alt_C.Cost = cell(nS,1);
    Alt_C.Binary = cell(nS,1);
    D.C = netDemand.C;
    for t = 1:1:nS %for every timestep
        [marginal,StorPower(t,:)] = instantMarginalCost(FirstProfile,scaleCost,netDemand,IC,dt,t);%update marginal cost
        QP = updateMatrices1Step(Plant.OneStep,Forecast,scaleCost(t,:),marginal,StorPower(t,:),dt,IC,FirstProfile,limit,t,T_build);
        QP_C = chillerOnly1step(QP,marginal,FirstProfile(t+1,:));
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Thermal Storage';}) && isfield(Plant.Generator(i).QPform.output,'C')
                D.C(t) = D.C(t) - StorPower(t,i);
            end
        end
        K_C = createCombinations(QP_C,D,FirstProfile(t+1,:),dt(t),t,[],[]);%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
        [BestDispatch,Alt_C] = eliminateCombinations(QP_C,K_C,Alt_C,IC,D,dt(t),t);%% combination elimination loop
        [GenOutput(t+1,:),IC] = updateIC(IC,BestDispatch,FirstProfile(t+1,:),dt(t),limit);
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
%     FirstProfile = checkAbChiller(FirstProfile,Forecast,dt);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'Thermal Storage';'Utility';'Electric Storage';'Solar';'Wind';})
            K_chill(:,i) = 1;  
        else
            K_chill(:,i) = FirstProfile(:,i)>0;
        end
    end
    IC = [CurrentState.Generators, CurrentState.Lines,CurrentState.Buildings(1,:)];
end

Prev = IC;
for t = 1:1:nS %for every timestep
    if ~isempty(FirstProfile) 
        [marginal,StorPower(t,:)] = instantMarginalCost(FirstProfile,scaleCost,netDemand,IC,dt,t);%update marginal cost
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
    else
        EC = [];
    end
    K_all = createCombinations(QP,netDemand,EC,dt(t),t,[],[]);%% create a matrix of all possible combinations 
    if isempty(K_chill)
        K = K_all;
    else
        K = createCombinations(QP,netDemand,EC,dt(t),t,K_chill(t+1,:),K_all);%% create a matrix of all possible combinations using the best chiller dispatch combination
    end
    [BestDispatch,Alt] = eliminateCombinations(QP,K,Alt,Prev,netDemand,dt(t),t);%% combination elimination loop
    j = 1;
    while ~isempty(Alt) && isempty(Alt.Disp{t}) && ~isempty(K_chill) && j<=length(Alt_C.Disp{t})
        K = createCombinations(QP,netDemand,EC,dt(t),t,Alt_C.Binary{t}(j,:),K_all);%% create a matrix of all possible combinations using the next best chiller dispatch combination
        [BestDispatch,Alt] = eliminateCombinations(QP,K,Alt,Prev,netDemand,dt(t),t);%% combination elimination loop
        j = j+1;
    end
    if ~isempty(FirstProfile) 
        [GenOutput(t+1,:),IC] = updateIC(IC,BestDispatch,FirstProfile(t+1,:),dt(t),limit);
    else
        [GenOutput,IC] = updateIC(IC,BestDispatch,[],dt(t),limit);
        GenOutput(1,1:nG) = GenOutput(1,1:nG) + Forecast.Renewable;
    end

     %update building temperatures
    if ~isempty(T_build)
        Tzone = zeros(2,nB);
        Twall = zeros(2,nB);
        for i = 1:1:nB
            [Tzone(:,i),Twall(:,i)] = BuildingSimulate(Plant.Building(i),Forecast.Weather.Tdb(t),Forecast.Weather.RH(t),dt(t)*3600,Forecast.Building.InternalGains(t,i),Forecast.Building.ExternalGains(t,i),Alt.Cooling{t}(i),Alt.Heating{t}(i),Forecast.Building.AirFlow(t,i),Forecast.Building.Damper(t,i),T_build(1,i),T_build(2,i));
        end
        T_build = [Tzone(2,:);Twall(2,:)];
    end
end
if ~isempty(FirstProfile)
    if isempty(K_chill)
        GenOutput = checkStartupCosts(GenOutput,StorPower,Alt,StartCost,dt,{'Electric Generator';'CHP Generator';'Heater';'Chiller';},40);
    else
        GenOutput = checkStartupCosts(GenOutput,StorPower,Alt,StartCost,dt,{'Electric Generator';'CHP Generator';'Heater';},40);
    end
    GenOutput(2:nS+1,1:nG) = GenOutput(2:nS+1,1:nG) + Forecast.Renewable;
end
% disp(['Time Spent in QP interations is ', num2str(sum(timeQP))]);
% disp(['Time not spent in QP is ', num2str(toc-sum(timeQP))]);
% disp(['Average # of combinations tested is ', num2str(mean(TestCombos))]);
end% Ends function StepByStepDispatch

function [GenOutput,IC] = updateIC(IC,BestDispatch,FirstProfile,dt,limit)
global Plant
nG = length(Plant.Generator);

if strcmp(limit,'unconstrained')
    IC = [];
    GenOutput = BestDispatch;
else
    EC = BestDispatch(1:nG);
    if ~isempty(FirstProfile)
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
                d_SOC = FirstProfile(i) - IC(i);%first profile already accounted for loss
                new_d_SOC = -BestDispatch(i)*dt/Plant.Generator(i).QPform.Stor.DischEff;
                EC(i) = IC(i) + d_SOC + new_d_SOC;
                if EC(i)<0
                    EC(i) = 0;
                    disp(strcat('Warning: ',Plant.Generator(i).Name,'_ is going negative'))
                elseif EC(i)>Plant.Generator(i).QPform.Stor.UsableSize
                    EC(i) = Plant.Generator(i).QPform.Stor.UsableSize;
                    disp(strcat('Warning: ',Plant.Generator(i).Name,'_ is exceeding max charge'))
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
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
                IC(i) = EC(i); %update the state of storage
            end
        end
    end
end
end%ends function updateIC