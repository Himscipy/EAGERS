function GenOutput = StepByStepDispatch(Forecast,scaleCost,dt,limit,FirstProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
global Plant CurrentState
nG = length(Plant.Generator);
nB = length(Plant.Building);
nL = length(Plant.OneStep.Organize.States)-nG-nB;
dX_dt = zeros(1,nG);
T = CurrentState.Buildings;
Dispatchable = logical(Plant.OneStep.Organize.Dispatchable);
if ~isempty(FirstProfile)
    IC = [CurrentState.Generators, CurrentState.Lines,CurrentState.Buildings];
    [nS,~] = size(scaleCost);
    StorPower = zeros(nS,nG);
    GenOutput = zeros(nS+1, nG+nL+nB);%nS should equal 1 in this case (finding IC)
    GenOutput(1,:) = IC;
    Alt.Disp = cell(nS,1);
    Alt.Cost = cell(nS,1);
    Alt.Binary = cell(nS,1);

    Binary = true(nS+1,nG);
    Binary(1,Dispatchable) = IC(Dispatchable)>0;
    uniqueCombinations = zeros(nS,1); %how many QP optimizations needed to be run in eliminate combinations
else
    nS = 1;
    StorPower = zeros(1,nG);
    IC = [];
end
UB = zeros(1,nG);
StartCost = zeros(1,nG);
Type = cell(nG,1);
QP.Organize.Enabled = zeros(1,nG);%message of which components are not enabled
nC = 0;
nH = 0;
for i = 1:1:nG
    if Plant.Generator(i).Enabled
        QP.Organize.Enabled(i) = 1;
        if Dispatchable(i) 
            Type(i) = {Plant.Generator(i).Type};
            if ismember(Plant.Generator(i).Type,{'Chiller'})
                nC = nC+1;
            elseif ismember(Plant.Generator(i).Type,{'Heater'})
                nH = nH+1;
            end
        end
    end
    if ~isempty(Plant.Generator(i).QPform.states)
        states = Plant.Generator(i).QPform.states(:,end);
        for j = 1:1:length(states)
            if strcmp((states{j}),'Y') && strcmp(Plant.Generator(i).Type,'Hydro Storage')
                %Need to convert UB to kW from kacre-ft ????
                UB(i) = UB(i) + ((12.1*Plant.Generator(i).QPform.(states{j}).ub(end))/Plant.Generator(i).QPform.Stor.Power2Flow); %kacr-ft to kw: 12.1*kacre-ft = kcfs..... kcfs/Power2Flow = kW
            elseif strcmp((states{j}),'S') && strcmp(Plant.Generator(i).Type,'Hydro Storage')
                %Don't add to UB; non-power producing
            else
                UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(end);
            end
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
        netDemand.H = 0;
        netDemand.C = 0;
        for i = 1:1:nB
            netDemand.E = netDemand.E + Forecast.Building(i).E0(t);
            netDemand.H = netDemand.H + Forecast.Building(i).H0(t);
            netDemand.C = netDemand.C + Forecast.Building(i).C0(t);
            Loads.Building(i).E0 = Forecast.Building(i).E0(t);
            Loads.Building(i).H0 = Forecast.Building(i).H0(t);
            Loads.Building(i).C0 = Forecast.Building(i).C0(t);
            Loads.Building(i).Tset_H = Forecast.Building(i).Tset_H(t);
            Loads.Building(i).Tset_C = Forecast.Building(i).Tset_C(t);
            Loads.Building(i).Tset = Forecast.Building(i).Tset(t);
        end
    end
    if isempty(FirstProfile) %finding initial conditions
        marginal = instantMarginalCost([],scaleCost);%update marginal cost
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable,scaleCost,marginal,dt,[],[],[],[],T);
    else
        if strcmp(limit, 'constrained')
            MinPower = max(0,IC(1:nG)-dX_dt*dt(t));
            MaxPower = min(UB,IC(1:nG)+dX_dt*dt(t));
        else
            MinPower = max(0,IC(1:nG)-dX_dt*sum(dt(1:t)));
            MaxPower = min(UB,IC(1:nG)+dX_dt*sum(dt(1:t)));
        end
        marginal = instantMarginalCost(FirstProfile(t+1,:),scaleCost(t,:));%update marginal cost
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
        QP = updateMatrices1Step(Plant.OneStep,Loads,Forecast.Renewable(t,:),scaleCost(t,:),marginal,dt(t),FirstProfile(t+1,:),StorPower(t,:),MinPower,MaxPower,T);
    end
    if nC>2 %more than 2 dispatchable chillers, perform seperate optimization and keep n_test best scenarios
        n_test = 3;%3 %number of chiller combinations to test when solving the generator problem
        D.C = netDemand.C;
        QP_C = chillerOnly1step(QP,marginal);
        K_C = createCombinations(QP_C,D,[],IC,dt(t));%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
        [ChillerDispatch,~,K_C] = eliminateCombinations(QP_C,K_C,n_test,D,dt(t));%% combination elimination loop
        netDemand = rmfield(netDemand,'C');
    else
        K_C = [];
    end
    %% Idea is to dispatch chillers with a electric and heating utility that has a price = marginal cost of generation, then find 3 best combinations
    %This can greatly reduce the combinations to be tested with the electric generators
    K = createCombinations(QP,netDemand,K_C,IC,dt(t));%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
    [lines,~] = size(K);
    
    [FeasibleDispatch,Cost,K] = eliminateCombinations(QP,K,lines,netDemand,dt(t));%% combination elimination loop
   
    if ~isempty(FirstProfile)
        Alt.Binary{t} = K>0;
        Alt.Disp{t} = FeasibleDispatch;
        uniqueCombinations(t) = length(Cost);
        if isempty(Alt.Disp{t})
            disp(['No feasible combination of generators at step' num2str(t)]);
            BestDispatch = IC;
            Binary(t+1,Dispatchable) = IC(Dispatchable)>0;
        else
            BestDispatch = Alt.Disp{t}(1,:);
            Alt.Cost{t} = Cost-Cost(1);
            Binary(t+1,:) = Alt.Binary{t}(1,:);
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
        if strcmp(limit, 'constrained')||strcmp(limit,'initially constrained')%if its constrained but not initially constrained then make the last output the initial condition
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
        GenOutput = FeasibleDispatch(1,:);
        GenOutput(:,1:nG) = GenOutput(:,1:nG) + Forecast.Renewable;
    end
end
if ~isempty(FirstProfile)
    GenOutput = checkStartupCosts(GenOutput,Alt,Binary,StartCost,dt,Type,{'Electric Generator';'CHP Generator';},20);
    if nC>0
        GenOutput = checkStartupCosts(GenOutput,Alt,Binary,StartCost,dt,Type,{'Chiller';},20);
    end
    GenOutput(2:nS+1,1:nG) = GenOutput(2:nS+1,1:nG) + Forecast.Renewable;
end
% disp(['Time Spent in QP interations is ', num2str(sum(timeQP))]);
% disp(['Time not spent in QP is ', num2str(toc-sum(timeQP))]);
% disp(['Average # of combinations tested is ', num2str(mean(TestCombos))]);
end% Ends function StepByStepDispatch