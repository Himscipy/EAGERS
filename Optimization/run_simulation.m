function [gen,buildings,cool_tower,design,dispatch,predicted] = run_simulation(date,num_steps,s_i,handles,test_data,hist_prof,gen,buildings,cool_tower,options,subnet,op_mat_a,op_mat_b,one_step,design,dispatch,predicted)
%num_steps: the number of dispatch optimizations that will occur during the entire simulation
%s_i: Counter for dispatch loop
freq = 1; %period of repetition (1 = 1 day)
res = options.Resolution/24;
n_o = round(freq/res)+1;
time = build_time_vector(options);%% set up vector of time interval
date = date+[0;time/24];
if num_steps == 0
    s_i = 1;
    [num_steps,dispatch,predicted,design,~] = pre_allocate_space(gen,buildings,cool_tower,subnet,options,test_data);
    % k = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
    k = 2;
    if k ==1
        [gen,buildings] = manual_ic(gen,buildings,cool_tower,subnet);
    else
        prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[]);
        now_data = get_data(test_data,date(1),[]);
        if strcmp(options.forecast,'Perfect')
            future_data = get_data(test_data,date(2:end),[]);
        else
            future_data = [];
        end
        [data_t0,gen,~] = update_forecast(gen,buildings,cool_tower,subnet,options,date(1),hist_prof,prev_data,now_data,future_data);
        [gen,cool_tower] = automatic_ic(gen,buildings,cool_tower,subnet,date(1),one_step,options,data_t0);% set the initial conditions
    end
    for i = 1:1:length(gen)
        dispatch.GeneratorState(1,i) = gen(i).CurrentState(1);
        if strcmp(gen(i).Type,'Hydro Storage')
            n = gen(i).QPform.Hydro.subnetNode;%dam #
            dispatch.hydroSOC(1,n) = gen(i).CurrentState(2);
        end
    end
    for i = 1:1:length(buildings)
        dispatch.Buildings(1,i) = buildings(i).Tzone;
    end
    dispatch.Timestamp(1) = date(1);
end
timers = zeros(num_steps,3); % To record times set to zeros(1,3), to not record set to [];
solution.Dispatch = [];
sim_waitbar=waitbar(0,'Running Dispatch','Visible','off');

while s_i<num_steps-1
    prev_data = get_data(test_data,linspace((date(1) - res - freq),date(1)-res,n_o)',[]);
    now_data = get_data(test_data,date(1),[]);
    if strcmp(options.forecast,'Perfect')
    	future_data = get_data(test_data,date(2:end),[]);
    else
        future_data = [];
    end
    [forecast,gen,buildings] = update_forecast(gen,buildings,cool_tower,subnet,options,date(2:end),hist_prof,prev_data,now_data,future_data);%% function that creates demand vector with time intervals coresponding to those selected
    solution = dispatch_loop(gen,buildings,cool_tower,subnet,op_mat_a,op_mat_b,one_step,options,date,forecast,solution);
    timers(s_i,:) = solution.timers;
    if ~isempty(handles)
        pause(0.01);
        stop = get(handles.Stop,'Value');
        if stop ==1
            return %stop button was pressed
        end
    end    
    
    [c,~,~] = net_cost(gen,solution.Dispatch,date,'Dispatch');
    predicted.Cost(s_i) = sum(c);
    [s_i,date,gen,buildings,cool_tower,design,dispatch,predicted] = dispatch_record(gen,buildings,cool_tower,subnet,options,test_data,s_i,date,forecast,solution,design,dispatch,predicted);
    if isfield(subnet,'Hydro')
        for n = 1:1:length(subnet.Hydro.nodes)
            test_data.Hydro.OutFlow(s_i,n) = dispatch.OutFlow(s_i,n);
        end
    end
    update_gui_status(handles,solution,date-options.Resolution/24,s_i-1,gen,options,dispatch)
    
    waitbar(s_i/num_steps,sim_waitbar,strcat('Running Dispatch'));
end