function solution = dispatch_loop(gen,building,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,last)
%% calculate optimal dispatch over the forecast horizon
dt = (date(2:end) - date(1:end-1))*24;
n_g= length(gen);
n_s = length(date)-1;

scale_cost = update_cost(date(2:end),gen); %% All feedstock costs were assumed to be 1 when building matrices 
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
predict_dispatch = zeros(length(date),n_g);
if isempty(last) || isempty(last.Dispatch)
    for i = 1:1:n_g
        predict_dispatch(:,i) = gen(i).CurrentState(1);
    end
else
    index = [(3:length(last.Dispatch(:,1)))';length(last.Dispatch(:,1));];
    for i = 1:1:n_g
        predict_dispatch(:,i) = [gen(i).CurrentState(1); last.Dispatch(index,i)];
    end
end
margin_cost = update_mc(gen,predict_dispatch,scale_cost,dt,[]);
if strcmp(options.solver,'Gurobi')
    qp = update_matrices(gen,building,cool_tower,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]);
    tic
    x = gurobi_opt(qp);
    solution = sort_solution(x,qp);
    tsim(1,1) = toc;
elseif strcmp(options.solver,'ANN')
    % Step 2: unit commitment, handled by a trained ANN
    tic
    %train the first time through
    training = isempty(last.Dispatch);
    locked = fireANN(last.Dispatch,forecast,scale_cost,training);
    qp = update_matrices(gen,building,cool_tower,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]);
    %make sure initial condition of locked is correct
    for i = 1:1:n_g
        if qp.Organize.Dispatchable(i) ==1
            locked(1,i) = gen(i).Status;
        end
    end
    tsim(1,2) = toc;
    % Step 3: complete optimization with fit B
    tic
    qp = disable_generators(qp,locked,[]);%Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution = sort_solution(x,qp);
    else %this means it was infeasible using the ANN, so use mcQP
        %default to mcQP if ANN doesn't work
        disp(strcat('ANN unit commitment infeasible, defaulting to mcQP at timestep', num2str(date(1))))
        solution = dispatch_loop(gen,building,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options.SpinReserve,true,'quadprog',date,forecast,last);
    end
    tsim(1,3) = toc;
else
    qp = update_matrices(gen,building,cool_tower,subnet,options,op_mat_a,date,scale_cost,margin_cost,forecast,[]);
    %% Step 1 Determine initial dispatch
    tic
    locked = true(n_s+1,n_g);
    for i = 1:1:n_g
        if ~gen(i).Enabled
            locked(:,i) = 0;
        end
    end
    qp = disable_generators(qp,locked,[]);%Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution1 = sort_solution(x,qp);
    else
        disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
        solution1.Dispatch = predict_dispatch;
    end
    tsim(1,1) = toc;
    if any(qp.Organize.Dispatchable) %might be some on/off combinations
        if options.MixedInteger
            tic
            %% Step 2:  dispatch step by step
            optimal_state = dispatch_step(gen,building,cool_tower,subnet,options,one_step,date(1),forecast,scale_cost,dt,solution1.Dispatch);
            clear mex
            tsim(1,2) = toc;
            tic
            %% Step 3:  2nd complete optimization
            margin_cost = update_mc(gen,optimal_state,scale_cost,dt,[]);
            qp_0 = update_matrices(gen,building,cool_tower,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]); %update fit B matrices
            locked = verify_ramping(gen,subnet,qp_0,locked,optimal_state,dt);
            qp = disable_generators(qp_0,locked,[]);%Disable generators here
            [x,feasible] = call_solver(qp);%this is the dispatch with fit B
            if feasible == 1
                solution = sort_solution(x,qp);
            else 
                disp('error: Cannot Find Feasible Dispatch');
                solution = solution1;
%                 [Solution, Feasible] = find_feasible(Gen,QP_0,Locked);
            end
        else
            tic
            margin_cost = update_mc(gen,solution1.Dispatch,scale_cost,dt,[]);
            qp_0 = update_matrices(gen,building,cool_tower,subnet,options,op_mat_b,date,scale_cost,margin_cost,forecast,[]); %update fit B matrices
            solution = cqp_method(gen,qp_0,solution1.Dispatch,date);
        end
        tsim(1,3) = toc;
    else
        solution = solution1;
    end 
end
solution.timers = tsim;
end%ends function dispatch_loop