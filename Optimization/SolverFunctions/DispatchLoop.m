function Solution = DispatchLoop(Date,Forecast,LastSolution)
%% calculate optimal dispatch over the forecast horizon
global Plant CurrentState %%  loaded by GUI & load generators
dt = (Date(2:end) - Date(1:end-1))*24;
nG = length(Plant.Generator);
nS = length(Date)-1;

scaleCost = updateGeneratorCost(Date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
if isempty(LastSolution) || isempty(LastSolution.Dispatch)
    PredictDispatch = ones(length(Date),1)*CurrentState.Generators;
else
    index = [(3:length(LastSolution.Dispatch(:,1)))';length(LastSolution.Dispatch(:,1));];
    PredictDispatch = [CurrentState.Generators; LastSolution.Dispatch(index,:)];
end
marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt,1);%the dispatch is whatever has been dispatched so far, except for the initial condition.

if strcmp(Plant.optimoptions.solver,'Gurobi')
    pause(10)
    tic
    QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
    x = gurobi_opt(QP);
    Solution = sortSolution(x,QP);
    tsim(1,1) = toc;
else
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
        Solution1 = sortSolution(x,QP);
    end
    tsim(1,1) = toc;
    if any(QP.Organize.Dispatchable) %might be some on/off combinations
        if ~(Feasible==1)%% hopefully not here
            disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
            Solution1.Dispatch = PredictDispatch;
        end
        %% Step 2:  dispatch step by step
        if Plant.optimoptions.MixedInteger
            tic
            OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,'initially constrained',Solution1.Dispatch);
            clear mex
            tsim(1,2) = toc;
        else
            OptimalState = Solution1.Dispatch;
        end
        %% Start with optimal dispatch, and check if feasible
        tic
        marginCost = updateMarginalCost(OptimalState,scaleCost,dt,2);
        QP_0 = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
        Locked = CheckRampRates(QP_0,Locked,OptimalState,dt);
        QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
        [x,Feasible] = callQPsolver(QP);%this is the dispatch with fit B
        if Feasible == 1
            Solution = sortSolution(x,QP);
        else 
            [Solution, Feasible] = FindFeasible(QP_0,Locked);
        end
        if Feasible==1
            if ~Plant.optimoptions.MixedInteger
                Solution = FilterGenerators(QP_0,Solution,Locked,Date);
                %add something to see if FilterGenerators changed anything
            end
        else
            disp('error: Cannot Find Feasible Dispatch');
        end
        tsim(1,3) = toc;
    end
end
Solution.timers = tsim;
end%End Dispatch Loop