function [SolutionGurobi,Solution_mcQP,Solution_cQP, fitADispatch, tsim] = GurobiTest(Date)
global Plant DateSim CurrentState
loadTestData
DateSim = Date;
Date = Date+[0;buildTimeVector(Plant.optimoptions)/24];
if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
    initializeOptimization
end
loadTestData
tsim = zeros(1,5);
reloadLast24hour(DateSim,Plant.optimoptions.Resolution)%re-load the previous 24 hours
interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Horizon/24,0.00);%create test data at correct frequency
Data = GetCurrentData(DateSim); 
automaticInitialCondition(Data);
Forecast = updateForecast(Date(2:end));
dt = (Date(2:end) - Date(1:end-1))*24;
nG = length(Plant.Generator);
nS = length(Date)-1;
scaleCost = updateGeneratorCost(Date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
PredictDispatch = ones(length(Date),1)*CurrentState.Generators;
marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt,[]);
Plant.optimoptions.MixedInteger = true;
QP = updateMatrices(Plant.OpMatA,Date,scaleCost,marginCost,Forecast,[]);
%% Step 1 Determine initial dispatch
Locked = true(nS+1,nG);
for i = 1:1:nG
    if ~Plant.Generator(i).Enabled
        Locked(:,i) = 0;
    end
end
tic
QP = disableGenerators(QP,Locked,[]);%Disable generators here
[x,Feasible] = callQPsolver(QP);
if Feasible == 1
    Solution = sortSolution(x,QP);
end
tsim(1,2) = toc;
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    Solution.Dispatch = PredictDispatch;
end
fitADispatch = Solution.Dispatch;

%% Step 2:  dispatch step by step
tic
OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,Solution.Dispatch);
clear mex
tsim(1,2) = toc;
%% Start with optimal dispatch, and check if feasible
tic
marginCost = updateMarginalCost(OptimalState,scaleCost,dt,[]);
QP_0 = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
Locked = CheckRampRates(QP_0,Locked,OptimalState,dt);
QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
[x,Feasible] = callQPsolver(QP);%this is the dispatch with fit B
if Feasible == 1
    Solution_mcQP = sortSolution(x,QP);
else
    Solution_mcQP = [];
end
tsim(1,4) = toc;


%% Run Gurobi
PredictDispatch = ones(length(Date),1)*CurrentState.Generators;
marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt,[]);
QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
QP.constCost = dt*QP.constCost.*scaleCost;
nG = length(Plant.Generator);

%add in end conditions so that storage is the same at the final timestep
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Stor')
        [m,n] = size(QP.Aeq);
        QP.Aeq = [QP.Aeq;zeros(1,n);];
        QP.Aeq(m+1,QP.organize{end,i}) = 1;
        QP.beq = [QP.beq;Solution_mcQP.Dispatch(end,i);];
    end
end
tic
x = gurobi_opt(QP);
if ~any(isnan(x))
    SolutionGurobi = sortSolution(x,QP);
else 
    SolutionGurobi = [];
end
tsim(1,1) = toc;


%% run cQP 
tic
[Solution_cQP,LBrelax] = cQP_Feasibility(Solution.Dispatch,Forecast,scaleCost,Date);
tsim(1,5) = toc;


%% compare costs
if ~isempty(Solution_cQP)
    [Cost,Eimbalance,Himbalance] = NetCostCalc(Solution_cQP.Dispatch,Date,'Dispatch');
    Solution_cQP.Cost = Cost;
    Solution_cQP.Eimbalance = Eimbalance;
    Solution_cQP.Himbalance = Himbalance;
    Solution_cQP.LBrelax = LBrelax;
end

if ~isempty(Solution_mcQP)
    [Cost,Eimbalance,Himbalance] = NetCostCalc(Solution_mcQP.Dispatch,Date,'Dispatch');
    Solution_mcQP.Cost = Cost;
    Solution_mcQP.Eimbalance = Eimbalance;
    Solution_mcQP.Himbalance = Himbalance;
end

if ~isempty(SolutionGurobi)
    [Cost,Eimbalance,Himbalance] = NetCostCalc(SolutionGurobi.Dispatch,Date,'Dispatch');
    SolutionGurobi.Cost = Cost;
    SolutionGurobi.Eimbalance = Eimbalance;
    SolutionGurobi.Himbalance = Himbalance;
end

