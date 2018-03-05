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
marginCost = updateMarginalCost(PredictDispatch,scaleCost,dt,[]);
if strcmp(Plant.optimoptions.solver,'Gurobi')
    QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
    tic
    x = gurobi_opt(QP);
    Solution = sortSolution(x,QP);
    tsim(1,1) = toc;
elseif strcmp(Plant.optimoptions.solver,'ANN')
    % Step 2: unit commitment, handled by a trained ANN
    tic
    %train the first time through
    training = isempty(LastSolution.Dispatch);
    Locked = fireANN(LastSolution.Dispatch,Forecast,scaleCost,training);
    QP = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]);
    %make sure initial condition of locked is correct
    for i = 1:1:nG
        if QP.Organize.Dispatchable(i) ==1
            Locked(1,i) = CurrentState.Generators(i)>0;
        end
    end
    tsim(1,2) = toc;
    % Step 3: complete optimization with fit B
    tic
    QP = disableGenerators(QP,Locked,[]);%Disable generators here
    [x,Feasible] = callQPsolver(QP);
    if Feasible == 1
        Solution = sortSolution(x,QP);
    else %this means it was infeasible using the ANN, so use mcQP
        %default to mcQP if ANN doesn't work
        disp(strcat('ANN unit commitment infeasible, defaulting to mcQP at timestep', num2str(Date(1))))
        Plant.optimoptions.MixedInteger = true;
        Plant.optimoptions.solver = 'quadprog';
        Solution = DispatchLoop(Date,Forecast,LastSolution);
        Plant.optimoptions.solver = 'ANN';
    end
    tsim(1,3) = toc;
else
    QP = updateMatrices(Plant.OpMatA,Date,scaleCost,marginCost,Forecast,[]);
    %% Step 1 Determine initial dispatch
    tic
    Locked = true(nS+1,nG);
    for i = 1:1:nG
        if ~Plant.Generator(i).Enabled
            Locked(:,i) = 0;
        end
    end
    QP = disableGenerators(QP,Locked,[]);%Disable generators here
    [x,Feasible] = callQPsolver(QP);
    if Feasible == 1
        Solution1 = sortSolution(x,QP);
    else
        disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
        Solution1.Dispatch = PredictDispatch;
    end
    tsim(1,1) = toc;
    if any(QP.Organize.Dispatchable) %might be some on/off combinations
        if Plant.optimoptions.MixedInteger
            tic
            %% Step 2:  dispatch step by step
            OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,Solution1.Dispatch);
            clear mex
            tsim(1,2) = toc;
            tic
            %% Step 3:  2nd complete optimization
            marginCost = updateMarginalCost(OptimalState,scaleCost,dt,[]);
            QP_0 = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
            Locked = CheckRampRates(QP_0,Locked,OptimalState,dt);
            QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
            [x,Feasible] = callQPsolver(QP);%this is the dispatch with fit B
            if Feasible == 1
                Solution = sortSolution(x,QP);
            else 
                [Solution, Feasible] = FindFeasible(QP_0,Locked);
            end
            if Feasible~=1
                disp('error: Cannot Find Feasible Dispatch');
                Solution = Solution1;
            end  
        else
            tic
            Solution = cQP_Feasibility(Solution1.Dispatch,Forecast,scaleCost,Date);
        end
        tsim(1,3) = toc;
    else
        Solution = Solution1;
    end 
end
Solution.timers = tsim;
end%End Dispatch Loop