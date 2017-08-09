function GenOutput = StepByStepDispatch(Forecast,scaleCost,dt,limit,FirstProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
global Plant CurrentState
if license('test','Distrib_Computing_Toolbox') 
    parallel = true;
else parallel = false;
end
nG = length(Plant.Generator);
nB = length(Plant.Building);
nL = length(Plant.OneStep.Organize.States)-nG-nB;
dX_dt = zeros(1,nG);
T = CurrentState.Buildings;
if ~isempty(FirstProfile)
    IC = [CurrentState.Generators, CurrentState.Lines,CurrentState.Buildings];
    [nS,~] = size(scaleCost);
    StorPower = zeros(nS,nG);
    GenOutput = zeros(nS+1, nG+nL+nB);%nS should equal 1 in this case (finding IC)
    GenOutput(1,:) = IC;
    I = zeros(nS,1);
    Alt.Disp = cell(nS,1);
    Alt.Cost = cell(nS,1);
    Alt.Binary = cell(nS,1);
    VentedHeat = zeros(nS,1);

    Binary = true(nS+1,nG);
    Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
    Binary(1,Dispatchable) = IC(Dispatchable)>0;
    uniqueCombinations = zeros(nS,1); %how many QP optimizations needed to be run in eliminate combinations
else
    nS = 1;
    StorPower = zeros(1,nG);
end
UB = zeros(1,nG);
StartCost = zeros(1,nG);
QP.Organize.Enabled = zeros(1,nG);%message of which components are not enabled
for i = 1:1:nG
    if Plant.Generator(i).Enabled
        QP.Organize.Enabled(i) = 1;
    end
    if~isempty(Plant.Generator(i).QPform.states)
        states = Plant.Generator(i).QPform.states(:,end);
        for j = 1:1:length(states);
            UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(end);
        end
    end
    if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
        StartCost(i) = Plant.Generator(i).VariableStruct.StartCost;
    end
    dX_dt(i) = Plant.Generator(i).VariableStruct.dX_dt;
end

if isfield(Forecast,'Demand')
    Outs = fieldnames(Forecast.Demand);
end
for t = 1:1:nS %for every timestep
    if isfield(Forecast,'Demand')
        for j = 1:1:length(Outs)
            Loads.(Outs{j}) = Forecast.Demand.(Outs{j})(t,:); %update demand
            netDemand.(Outs{j}) = sum(Loads.(Outs{j}));
        end
    else 
        netDemand.E = 0;
        for i = 1:1:nB
            netDemand.E = netDemand.E + Forecast.Building(i).E0(t);
            Loads.Building(i).E0 = Forecast.Building(i).E0(t);
            Loads.Building(i).H0 = Forecast.Building(i).H0(t);
            Loads.Building(i).C0 = Forecast.Building(i).C0(t);
            Loads.Building(i).Tset_H = Forecast.Building(i).Tset_H(t);
            Loads.Building(i).Tset_C = Forecast.Building(i).Tset_C(t);
            Loads.Building(i).Tset = Forecast.Building(i).Tset(t);
        end
    end
    if isempty(FirstProfile) %finding initial conditions
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable,scaleCost,dt,[],[],[],[],T);
    else
        if strcmp(limit, 'constrained')
            MinPower = max(0,IC(1:nG)-dX_dt*dt(t));
            MaxPower = min(UB,IC(1:nG)+dX_dt*dt(t));
        else
            MinPower = max(0,IC(1:nG)-dX_dt*sum(dt(1:t)));
            MaxPower = min(UB,IC(1:nG)+dX_dt*sum(dt(1:t)));
        end
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';}) 
                loss = dt(t)*(Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
                Power = (IC(i) - FirstProfile(t+1,i) + loss)/dt(t);%expected output of storage in kW to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
                if Power>0 %discharging
                    StorPower(t,i) = Power*Plant.Generator(i).QPform.Stor.DischEff;
                else %charging
                    StorPower(t,i) = Power/Plant.Generator(i).QPform.Stor.ChargeEff; 
                end
            elseif ismember(Plant.Generator(i).Type,{'Hydro Storage';}) 
                %lineout is the downriver segment associated with this dam
                StorPower(t,i) = FirstProfile(t+1,i);
                n = Plant.Generator(i).QPform.subnetNode;
                lineOut = Plant.subNet.Hydro.lineNumber(n);
                if t>1
                    FirstProfile(t+1,lineOut) = FirstProfile(t+1,lineOut) - BestDispatch(lineOut); %subtract excess flow from last timestep
                end
            end
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
                out = fieldnames(Plant.Generator(i).QPform.output);
                netDemand.(out{1}) = netDemand.(out{1}) - StorPower(t,i);
            end
        end
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable(t,:),scaleCost(t,:),dt(t),FirstProfile(t+1,:),StorPower(t,:),MinPower,MaxPower,T);
    end

    K = createCombinations(QP,netDemand);%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
    [lines,~] = size(K);
    if lines == 0
        disp('Zero feasible combinations to test: ERROR')
    end
    FeasibleDispatch = zeros(lines,nG+nL+nB);
    HeatVent = zeros(lines,1);%net heat loss
    Cost = zeros(lines,1);
    feasible = false(lines,1);
    if parallel
        parfor i = 1:lines
            [FeasibleDispatch(i,:),Cost(i),feasible(i),HeatVent(i),~] = eliminateCombinations(QP,netDemand,K(i,:),parallel,[],[],[]);%% combination elimination loop
        end
    else
        nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
        [~, line] = sort(nzK); %sort the rows by number of generators that are on
        K = K(line,:);
        %% test the cases for cost
        for i=1:lines %run the quadprog/linprog for all cases with the least number of generators
            if i<=length(K(:,1))
                [FeasibleDispatch(i,:),Cost(i),feasible(i),HeatVent(i),K] = eliminateCombinations(QP,netDemand,K(i,:),parallel,K,i,min(Cost(1:i)),dt(t));%% combination elimination loop
            end
        end
    end
    Cost = Cost(feasible);
    if isempty(Cost)
        disp('Zero feasible outcomes in StepByStep: ERROR')
    end
    [~,I(t)] = min(Cost);
    if ~isempty(FirstProfile)
        Alt.Binary{t} = K(feasible,:)>0;
        Alt.Disp{t} = FeasibleDispatch(feasible,:);
        uniqueCombinations(t) = length(Cost);
        VentedHeat(t) = HeatVent(I(t));
        if isempty(Alt.Disp{t})
            disp(['No feasible combination of generators at step' num2str(t)]);
            BestDispatch = IC;
            Binary(t+1,Dispatchable) = IC(Dispatchable)>0;
        else
            BestDispatch = Alt.Disp{t}(I(t),:);
            Alt.Cost{t} = Cost-Cost(I(t));
            Binary(t+1,:) = Alt.Binary{t}(I(t),:);
        end
        EC = BestDispatch(1:nG);
        for i = 1:1:nG
            if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
                Energy = (BestDispatch(i)+StorPower(t,i))*dt(t);
                loss = dt(t)*(Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
                if Energy>0 %discharging
                    EC(i) = IC(i) - Energy/Plant.Generator(i).QPform.Stor.DischEff - loss;%change in storage for this power output
                else %charging
                    EC(i) = IC(i) - Energy*Plant.Generator(i).QPform.Stor.ChargeEff - loss;%change in storage for this power output
                end
            end
        end
        if strcmp(limit, 'constrained')%if its constrained but not initially constrained then make the last output the initial condition
            IC = EC;
        else
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'})
                    IC(i) = EC(i); %update the state of storage
                end
            end
        end
        %update building temperatures
        if ~isempty(T)
            %% Do something here
        end
        GenOutput(t+1,1:nG) = EC;
        if nL>0||nB>0
            GenOutput(t+1,nG+1:nG+nL+nB) = BestDispatch(nG+1:nG+nL+nB);
        end
    else
        feas = nonzeros(feasible.*(1:1:lines)');
        GenOutput = FeasibleDispatch(feas(I),:);
        GenOutput(:,1:nG) = GenOutput(:,1:nG) + Forecast.Renewable;
    end
end
if ~isempty(FirstProfile)
    GenOutput(2:nS+1,:) = checkStartupCosts(Alt,Binary,StartCost,dt,I,nL+nB);
    GenOutput(2:nS+1,:) = GenOutput(2:nS+1,:) + Forecast.Renewable;
end
% disp(['Time Spent in QP interations is ', num2str(sum(timeQP))]);
% disp(['Time not spent in QP is ', num2str(toc-sum(timeQP))]);
% disp(['Average # of combinations tested is ', num2str(mean(TestCombos))]);
end% Ends function StepByStepDispatch