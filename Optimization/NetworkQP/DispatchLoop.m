function [GenDisp,tsim] = DispatchLoop(Date,Forecast,LastDispatch)
%% calculate optimal dispatch over the forecast horizon
global Plant CurrentState %%  loaded by GUI & load generators
dt = (Date(2:end) - Date(1:end-1))*24;
nG = length(Plant.Generator);
nS = length(Date)-1;

%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
if isempty(LastDispatch)
    PredictDispatch = ones(length(Date),1)*[CurrentState.Generators, CurrentState.Lines, CurrentState.Buildings];
else
    PredictDispatch = [CurrentState.Generators, CurrentState.Lines, CurrentState.Buildings; LastDispatch(3:end,:);LastDispatch(end,:)];
end
scaleCost = updateGeneratorCost(Date(2:end)); %% All feedstock costs were assumed to be 1 when building matrices 
if strcmp(Plant.optimoptions.solver,'NREL')
    tic
    GenDisp = NRELoptimization(Forecast,scaleCost);
    tsim(1,1) = toc;
else
    marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt);%the dispatch is whatever has been dispatched so far, except for the initial condition.
    Plant.OpMatA.Renewable = Forecast.Renewable;
    QP = updateMatrices(Plant.OpMatA,Date,scaleCost,marginCost,Forecast,[]);

    %% Step 1 Determine initial dispatch
    Locked = true(nS+1,nG);
    for i = 1:1:nG
        if ~Plant.Generator(i).Enabled
            Locked(:,i) = 0;
        end
    end
    tic
    [FirstDisp,Feasible] = DispatchQP(QP,Locked);
    tsim(1,1) = toc;

    if ~(Feasible==1)%% hopefully not here
        disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
        FirstDisp= PredictDispatch;
    end
    %% Step 2:  dispatch step by step
    if Plant.optimoptions.MixedInteger
        tic
        OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,'initially constrained',FirstDisp);
        clear mex
        tsim(1,2) = toc;
    else
        OptimalState = FirstDisp;
    end
    %% Start with optimal dispatch, and check if feasible
    tic
    marginCost = updateMarginalCost(OptimalState,scaleCost,dt);
    Plant.OpMatB.Renewable = Forecast.Renewable;
    QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
    for i = 1:1:nG
        if QP.Organize.Dispatchable(i) ==1
            Locked(OptimalState(:,i)==0,i)=false;
        end
    end
    [GenDisp, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B
    if Feasible ~=1
        [GenDisp, Feasible] = FindFeasible(QP,Locked);
    end
    if Feasible==1
        if ~Plant.optimoptions.MixedInteger
            GenDisp = FilterGenerators(QP,GenDisp,Locked,Date);
            %add something to see if FilterGenerators changed anything
        end
    else
        disp('error: Cannot Find Feasible Dispatch');
    end
    tsim(1,3) = toc;
    % Cost = sum(NetCostCalc(GenDisp,[0;Time/24]+Date,'Dispatch'));
end