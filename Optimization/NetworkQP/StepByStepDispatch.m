function GenOutput = StepByStepDispatch(Forecast,scaleCost,dt,IC,limit,FirstProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
global Plant
if license('test','Distrib_Computing_Toolbox') 
    parallel = true;
else parallel = false;
end
nG = length(Plant.Generator);
dX_dt = zeros(1,nG);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG;
ic = (~isempty(IC));
UB = zeros(1,nG);
StartCost = zeros(1,nG);
stor = [];
for i = 1:1:nG
    states = Plant.Generator(i).OpMatB.states;
    for j = 1:1:length(states);
        UB(i) = UB(i) + Plant.Generator(i).OpMatB.(states{j}).ub;
    end
    if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
        StartCost(i) = Plant.Generator(i).VariableStruct.StartCost;
    end
    if strfind(Plant.Generator(i).Type, 'Storage')
        stor(end+1) = i;
    end
    dX_dt(i) = Plant.Generator(i).VariableStruct.dX_dt;
end
if isempty(IC)
    nS = 0;
    StorPower = zeros(1,nG);
else
    [nS,~] = size(scaleCost);
    StorPower = zeros(nS,nG);
end
GenOutput = zeros(nS+1, nG+nL);%nS should equal 1 in this case (finding IC)
if ~isempty(IC)
    GenOutput(1,:) = IC;
end
TestCombos = zeros(nS,1); %how many QP optimizations needed to be run in eliminate combinations
Iterate = zeros(nS,1);
timeQP = zeros(nS,1);
%% note #2:  need to account for self discharging of storage
Outs = fieldnames(Forecast.Demand);
% Outs = Outs(~strcmp('T',Outs));
% Outs = Outs(~strcmp('Hydro',Outs));
I = zeros(nS,1);
Alt.Disp = cell(nS,1);
Alt.Cost = cell(nS,1);
Alt.Binary = cell(nS,1);
VentedHeat = zeros(nS,1);
Binary = true(nS+1,nG);
Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
if ic
    Binary(1,Dispatchable) = IC(Dispatchable)>0;
end
for t = 1:1:max(1,nS) %for every timestep
    for j = 1:1:length(Outs)
        Loads.(Outs{j}) = Forecast.Demand.(Outs{j})(t,:); %update demand
        netDemand.(Outs{j}) = sum(Loads.(Outs{j}));
    end
    if isfield(Plant.subNet,'Hydro')
        for i = 1:1:length(Plant.subNet.Hydro) %run through all the nodes in this network
            %% Node inflows (source/sink terms and upstream flow at time t-T ago if t<T)
            %need to be able to forecast sink/source
            I = find(strcmp(Plant.subNet.Hydro(i).nodes{1},Plant.Data.Hydro.Nodes));%column index of this node in the stored matrices of Forecast.Hydro.SourceSink, Forecast.Hydro.SpillFlow and Forecast.Hydro.Inflow
            equip = Plant.subNet.Hydro(i).Equipment;
            for k = 1:1:length(equip)
                if strcmp(Plant.Generator(equip(k)).Type,'Hydro Storage')
                    if isempty(FirstProfile) %finding initial conditions
                        StorPower(t,equip(k)) = (Forecast.Hydro.OutFlow(I,1) - Forecast.Hydro.SpillFlow(I,1))*Plant.Generator(equip(k)).OpMatA.output.E;
                    else
                        FlowGen = (IC(equip(k)) - FirstProfile(t+1,equip(k)))*12.1/dt(t) + Forecast.Hydro.SourceSink(I,t) + Forecast.Hydro.InFlow(I,t) - Forecast.Hydro.SpillFlow(I,t); %flow through turbines in 100 cfs
                        StorPower(t,equip(k)) = FlowGen*Plant.Generator(equip(k)).OpMatA.output.E;
                    end
                end
            end
        end
    end
    if isempty(FirstProfile) %finding initial conditions
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable,scaleCost(t,:),dt(t),[],StorPower,[],[]);
    else
        if strcmp(limit, 'constrained')
            MinPower = max(0,IC(1:nG)-dX_dt*dt(t));
            MaxPower = min(UB,IC(1:nG)+dX_dt*dt(t));
        else
            MinPower = max(0,IC(1:nG)-dX_dt*sum(dt(1:t)));
            MaxPower = min(UB,IC(1:nG)+dX_dt*sum(dt(1:t)));
        end
        for i = 1:1:length(stor)
            if ~strcmp(Plant.Generator(stor(i)).Type,'Hydro Storage') %hydro is organized by node, not generator
                loss = dt(t)*(Plant.Generator(stor(i)).OpMatA.Stor.SelfDischarge*Plant.Generator(stor(i)).OpMatA.Stor.UsableSize);
                Power = (IC(stor(i)) - FirstProfile(t+1,stor(i)) + loss)/dt(t);%expected output of storage in kW to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
                if Power>0 %discharging
                    StorPower(t,stor(i)) = Power*Plant.Generator(stor(i)).OpMatA.Stor.DischEff;
                else %charging
                    StorPower(t,stor(i)) = Power/Plant.Generator(stor(i)).OpMatA.Stor.ChargeEff; 
                end
            end
        end
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable(t,:),scaleCost(t,:),dt(t),FirstProfile(t,:),StorPower(t,:),MinPower,MaxPower);
    end
    QP.Organize.Enabled = zeros(1,nG);%message of which components are not enabled
    for i = 1:1:nG
        if Plant.Generator(i).Enabled
            QP.Organize.Enabled(i) = 1;
        end
    end
    timeQP(t) = toc;
    K = createCombinations(QP,netDemand);%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
    [lines,~] = size(K);
    if lines == 0
        disp('Zero feasible combinations to test: ERROR')
    end
    FeasibleDispatch = zeros(lines,nG+nL);
    HeatVent = zeros(lines,1);
    Cost = zeros(lines,1);
    feasible = false(lines,1);
    Iterations = zeros(lines,1);
    if parallel
        parfor i = 1:lines
            [FeasibleDispatch(i,:),Cost(i),feasible(i),Iterations(i),HeatVent(i),~] = eliminateCombinations(QP,netDemand,K(i,:),parallel,[],[],[]);%% combination elimination loop
        end
    else
        nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
        [~, line] = sort(nzK); %sort the rows by number of generators that are on
        K = K(line,:);
        %% test the cases for cost
        for i=1:lines %run the quadprog/linprog for all cases with the least number of generators
            if i<=length(K(:,1))
                [FeasibleDispatch(i,:),Cost(i),feasible(i),Iterations(i),HeatVent(i),K] = eliminateCombinations(QP,netDemand,K(i,:),parallel,K,i,min(Cost(1:i)),dt(t));%% combination elimination loop
            end
        end
    end
    timeQP(t) = toc - timeQP(t);
    Cost = Cost(feasible);
    if isempty(Cost)
        disp('Zero feasible outcomes in StepByStep: ERROR')
    end
    [~,I(t)] = min(Cost);
    Alt.Binary{t} = K(feasible,:)>0;
    Alt.Disp{t} = FeasibleDispatch(feasible,:);
    TestCombos(t) = length(Cost);
    Iterate(t) = sum(Iterations(feasible))/TestCombos(t);
    VentedHeat(t) = HeatVent(I(t));
    if isempty(Alt.Disp{t})
        disp(['No feasible combination of generators at step' num2str(t)]);
        BestDispatch = [IC,zeros(1,nL)];
        BestDispatch(stor) = 0;
        Binary(t+ic,Dispatchable) = IC(Dispatchable)>0;
    else
        BestDispatch = Alt.Disp{t}(I(t),:);
        Alt.Cost{t} = Cost-Cost(I(t));
        Binary(t+ic,:) = Alt.Binary{t}(I(t),:);
    end
    EC = BestDispatch(1:nG);
    if ~isempty(IC)
        for i = 1:1:length(stor)
            loss = dt(t)*(Plant.Generator(stor(i)).OpMatA.Stor.SelfDischarge*Plant.Generator(stor(i)).OpMatA.Stor.UsableSize);
            if strcmp(Plant.Generator(stor(i)).Type,'Hydro Storage')
                FlowGen = (BestDispatch(stor(i))+StorPower(t,stor(i)))/Plant.Generator(stor(i)).OpMatA.output.E; %conversion from power (kw) to mass flow (1000 ft^3/s)
                EC(stor(i)) = IC(stor(i)) + (Forecast.Hydro.InFlow(I,t) + Forecast.Hydro.SourceSink(I,t) - FlowGen - Forecast.Hydro.SpillFlow(I,t))*dt(t)/12.1;
            else
                Energy = (BestDispatch(stor(i))+StorPower(t,stor(i)))*dt(t);
                if Energy>0 %discharging
                    EC(stor(i)) = IC(stor(i)) - Energy/Plant.Generator(stor(i)).OpMatA.Stor.DischEff - loss;%change in storage for this power output
                else %charging
                    EC(stor(i)) = IC(stor(i)) - Energy*Plant.Generator(stor(i)).OpMatA.Stor.ChargeEff - loss;%change in storage for this power output
                end
            end
        end
        if strcmp(limit, 'constrained')%if its constrained but not initially constrained then make the last output the initial condition
            IC = EC;
        else
            IC(stor) = EC(stor); %update the state of storage
        end
    end
    GenOutput(t+(~isempty(IC)),1:nG) = EC;
    if nL>0
        GenOutput(t+(~isempty(IC)),nG+1:nG+nL) = BestDispatch(nG+1:nG+nL);
    end
end
if ic 
    GenOutput(2:nS+1,:) = checkStartupCosts(Alt,Binary,StartCost,dt,I,nL);
end
for i = 1:1:nG
    if any(Forecast.Renewable(:,i)>0)
        GenOutput(1+ic:end,i) = Forecast.Renewable(:,i);
    end
end
% disp(['Time Spent in QP interations is ', num2str(sum(timeQP))]);
% disp(['Time not spent in QP is ', num2str(toc-sum(timeQP))]);
% disp(['Average iterations for each QP is ', num2str(mean(Iterate))]);
% disp(['Average # of combinations tested is ', num2str(mean(TestCombos))]);