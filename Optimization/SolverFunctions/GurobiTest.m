function [SolutionGurobi,Solution_QP] = GurobiTest(Date)
global Plant DateSim CurrentState RealTimeData
DateSim = Date;
Date = Date+[0;buildTimeVector(Plant.optimoptions)/24];
if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
    initializeOptimization
end
loadTestData
resetLast24Hour
timestamp = (DateSim:Plant.optimoptions.Resolution/24:(DateSim+Plant.optimoptions.Horizon/24))';
RealTimeData = interpolateData(timestamp,0.00);%create test data at correct frequency
if isempty(Plant.Building)
    Data = GetHistoricalData(DateSim); 
else
    Data = updateForecast(DateSim,[]);
end
automaticInitialCondition(Data);
Forecast = updateForecast(Date(2:end),Data);
dt = (Date(2:end) - Date(1:end-1))*24;
nG = length(Plant.Generator);
nS = length(Date)-1;
scaleCost = updateGeneratorCost(Date(2:end)); %% All feedstock costs were assumed to be 1 when building matrices 
PredictDispatch = ones(length(Date),1)*[CurrentState.Generators, CurrentState.Lines, CurrentState.Buildings];
marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt,1);%the dispatch is whatever has been dispatched so far, except for the initial condition.
tic
QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
x = gurobi_opt(QP,dt);
if ~any(isnan(x))
    SolutionGurobi = sortSolution(x,QP);
else 
    SolutionGurobi = [];
end
tsim(1,1) = toc;

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
%% Step 2:  dispatch step by step
tic
OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,'initially constrained',Solution.Dispatch);
clear mex
tsim(1,3) = toc;
%% Start with optimal dispatch, and check if feasible
tic
marginCost = updateMarginalCost(OptimalState,scaleCost,dt,2);
QP_0 = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
Locked = CheckRampRates(QP_0,Locked,OptimalState,dt);
QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
[x,Feasible] = callQPsolver(QP);%this is the dispatch with fit B
if Feasible == 1
    Solution_QP = sortSolution(x,QP);
else
    Solution_QP = [];
end
tsim(1,4) = toc;