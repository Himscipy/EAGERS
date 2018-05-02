function gen_output = dispatch_step(gen,building,cool_tower,subnet,options,one_step,date,forecast,scale_cost,dt,first_profile)
% time is the time from the current simulation time, positive numbers are forward looking
% stor_power is the amount of power coming from (positive) or going into (negative) the storage device at each timestep according to the first dispatch
n_g = length(gen);
n_b = length(building);
n_ct = length(cool_tower);
temperatures.build = zeros(2,n_b);
for i = 1:1:n_b
    temperatures.build(1,i) = building(i).Tzone;
    temperatures.build(2,i) = building(i).Twall;
end
temperatures.cool_tower = zeros(1,n_ct);
for i = 1:1:n_ct
    temperatures.cool_tower(i) = cool_tower(i).fluid_temperature;
end
net_demand = agregate_demand(forecast);
limit = 'initially constrained';
start_cost = zeros(1,n_g);
ic = zeros(1,n_g);
for i = 1:1:n_g
    if isfield(gen(i).VariableStruct, 'StartCost')
        start_cost(i) = gen(i).VariableStruct.StartCost;
    end
    ic(i) = gen(i).CurrentState;
end
[nS,~] = size(scale_cost);
gen_output = zeros(nS+1,length(ic));
gen_output(1,:) = ic;
alt.Disp = cell(nS,1);
alt.Cost = cell(nS,1);
alt.Binary = cell(nS,1);

%If there are more than 2 chillers, solve it as a seperate problem
[combinations_chill,alt_c,first_profile] = seperate_chiller_problem(gen,building,cool_tower,subnet,options,one_step,date,forecast,first_profile,net_demand,scale_cost,start_cost,ic,dt,temperatures);
stor_power = zeros(nS,n_g);
for t = 1:1:nS %for every timestep
    if isfield(net_demand,'H')
        v_h = value_heat(gen,[ic;first_profile(t+1,:)],net_demand.H(t),dt(t));
    else
        v_h = value_heat(gen,[ic;first_profile(t+1,:)],0,dt(t));
    end
    stor_power(t,:) = find_stor_power(gen,[ic;first_profile(t+1,:)],dt(t));
    marginal = update_mc(gen,first_profile(t+1,:),scale_cost(t,:),dt(t),v_h);%update marginal cost
    QP = update_matrices_step(gen,building,cool_tower,subnet,options,one_step,forecast,scale_cost(t,:),marginal,stor_power(t,:),dt,ic,first_profile,limit,t,temperatures);
    QP.Organize.Enabled = true(1,n_g);%which components are enabled
    for i = 1:1:n_g
        if ismember(gen(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
            out = fieldnames(gen(i).QPform.output);
            if isfield(net_demand,out{1})
                net_demand.(out{1})(t) = net_demand.(out{1})(t) - stor_power(t,i);
            end
        end
        if ~gen(i).Enabled
            QP.Organize.Enabled(i) = false;
        end
    end
    ec = first_profile(t+1,:);
    prev = gen_output(t,:);
    combinations = create_combinations(gen,options,QP,net_demand,ec,dt(t),t,combinations_chill(t+1,:),'E');%% create a matrix of all possible combinations using the best chiller dispatch combination
    [best_dispatch,alt] = check_combinations(QP,combinations,alt,prev,net_demand,dt(t),t);%% combination elimination loop

    %If there was not actually a feasible combination of generators with the best chiller combination, try the next best chiller combination
    j = 1;
    while isempty(alt.Disp{t}) && ~isempty(alt_c) && j<=length(alt_c.Disp{t}(:,1))
        combinations = create_combinations(gen,options,QP,net_demand,ec,dt(t),t,alt_c.Binary{t}(j,:),'E');%% create a matrix of all possible combinations using the next best chiller dispatch combination
        [best_dispatch,alt] = check_combinations(QP,combinations,alt,prev,net_demand,dt(t),t);%% combination elimination loop
        j = j+1;
    end

    %update Initial conditions and building temperatures
    [gen_output(t+1,:),ic] = update_ic(gen,ic,best_dispatch,first_profile(t+1,:),dt(t),limit);%only updates the storage states in IC
    for i = 1:1:n_b
        [zone,wall] = building_simulate(building(i),...
            forecast.Weather.Tdb(t),forecast.Weather.RH(t),dt(t)*3600,...
            forecast.Building.InternalGains(t,i),forecast.Building.ExternalGains(t,i),...
            alt.Cooling{t}(i),alt.Heating{t}(i),forecast.Building.AirFlow(t,i),...
            forecast.Building.Damper(t,i),temperatures.build(1,i),temperatures.build(2,i));
        temperatures.build(:,i) = [zone(2);wall(2)];
    end
    for i = 1:1:n_ct
        imbalance = 0;
        equip = subnet.CoolingWater.Equipment{i};
        for k = 1:1:length(equip)
            j = equip(k);
            if strcmp(gen(j).Type,'Chiller')
                imbalance = imbalance + best_dispatch(j) + chill_input(gen(j).QPform,best_dispatch(j));
            elseif strcmp(gen(j).Type,'Cooling Tower')
                imbalance = imbalance - best_dispatch(j);
            end
        end
        capacitance = cool_tower(i).fluid_capacity*cool_tower(i).fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
        temperatures.cool_tower(i) = temperatures.cool_tower(i) + dt(t)*3600/capacitance*imbalance;
    end
end

%change best dispatch based on re-start costs
include = {'Electric Generator';'CHP Generator';'Heater';'Electrolyzer';'Hydrogen Generator';'Cooling Tower';};
if isempty(alt_c)
    include(end+1) = {'Chiller'};
end
dispatchable = logical(one_step.Organize.Dispatchable);
gen_output = minimize_start_costs(gen,[date;forecast.Timestamp],dispatchable,gen_output,stor_power,alt,start_cost,dt,include,40);
if isfield(forecast,'Renewable')
    gen_output(2:nS+1,1:n_g) = gen_output(2:nS+1,1:n_g) + forecast.Renewable;
end
end% Ends function dispatch_step

function [gen_output,ic] = update_ic(gen,ic,best_dispatch,first_profile,dt,limit)
n_g = length(gen);
gen_output = best_dispatch;
if strcmp(limit,'unconstrained')
    ic = [];
else
    ec = best_dispatch(1:n_g);
    if ~isempty(first_profile)
        for i = 1:1:n_g
            if ismember(gen(i).Type,{'Electric Storage';'Thermal Storage'})
                d_SOC = first_profile(i) - ic(i);%first profile already accounted for loss
                new_d_SOC = -best_dispatch(i)*dt/gen(i).QPform.Stor.DischEff;
                ec(i) = ic(i) + d_SOC + new_d_SOC;
                if ec(i)<0
                    error = ec(i)/gen(i).QPform.Stor.UsableSize*100;
                    ec(i) = 0;
                    disp(strcat('Warning: ',gen(i).Name,'_ is going negative by_',num2str(error),'%'))
                elseif ec(i)>gen(i).QPform.Stor.UsableSize
                    error = (ec(i)-gen(i).QPform.Stor.UsableSize)/gen(i).QPform.Stor.UsableSize*100;
                    ec(i) = gen(i).QPform.Stor.UsableSize;
                    disp(strcat('Warning: ',gen(i).Name,'_ is exceeding max charge by_',num2str(error),'%'))
                end
            end
        end
    end
    gen_output(1,1:n_g) = ec;
    gen_output(1,n_g+1:end) = best_dispatch(n_g+1:end);
    if strcmp(limit,'constrained')%if its constrained but not initially constrained then make the last output the initial condition
        ic = ec;
    else %if strcmp(limit,'initially constrained')
        for i = 1:1:n_g
            if isfield(gen(i).QPform,'Stor')
                ic(i) = ec(i); %update the state of storage
            end
        end
    end
end
end%ends function update_ic

function [combinations_chill,alt_c,first_profile] = seperate_chiller_problem(gen,building,cool_tower,subnet,options,one_step,date,forecast,first_profile,net_demand,scale_cost,start_cost,ic,dt,temperatures)
%%optimize chilling dispatch first
alt_c = [];
n_c = 0;
dispatchable = logical(one_step.Organize.Dispatchable);
n_g = length(gen);
for i = 1:1:n_g
    if gen(i).Enabled && dispatchable(i) 
        if ismember(gen(i).Type,{'Chiller'})
            n_c = n_c+1;
        end
    end
end
[n_s,~] = size(scale_cost);
combinations_chill = zeros(n_s+1,0);
stor_power = zeros(n_s,n_g);
if n_c>2 && any(net_demand.C>0)%more than 2 dispatchable chillers, perform seperate optimization and keep n_test best scenarios
    gen_output = zeros(n_s+1,length(ic));
    gen_output(1,:) = ic;
    combinations_chill = zeros(n_s+1,n_g);
    alt_c.Disp = cell(n_s,1);
    alt_c.Cost = cell(n_s,1);
    alt_c.Binary = cell(n_s,1);
    dem.C = net_demand.C;
    for t = 1:1:n_s %for every timestep
        if isfield(net_demand,'H')
            v_h = value_heat(gen,[ic;first_profile(t+1,:)],net_demand.H(t),dt(t));
        else
            v_h = value_heat(gen,[ic;first_profile(t+1,:)],0,dt(t));
        end
        stor_power(t,:) = find_stor_power(gen,[ic;first_profile(t+1,:)],dt(t));
        marginal = update_mc(gen,first_profile(t+1,:),scale_cost(t,:),dt(t),v_h);%update marginal cost
        QP = update_matrices_step(gen,building,cool_tower,subnet,options,one_step,forecast,scale_cost(t,:),marginal,stor_power(t,:),dt,ic,first_profile,'initially constrained',t,temperatures);
        QP_C = chiller_step(gen,building,cool_tower,subnet,options.SpinReserve,QP,marginal,first_profile(t+1,:));
        QP_C.Organize.Enabled = true(1,n_g);%which components are enabled
        for i = 1:1:n_g
            if ismember(gen(i).Type,{'Thermal Storage';}) && isfield(gen(i).QPform.output,'C')
                dem.C(t) = dem.C(t) - stor_power(t,i);
            end
            if ~gen(i).Enabled
                QP_C.Organize.Enabled(i) = false;
            end
        end
        combinations = create_combinations(gen,options,QP_C,dem,first_profile(t+1,:),dt(t),t,[],'C');%% create a matrix of all possible combinations (keep electrical and heating together if there are CHP generators, otherwise seperate by product)
        [best_dispatch,alt_c] = check_combinations(QP_C,combinations,alt_c,ic,dem,dt(t),t);%% combination elimination loop
        [gen_output(t+1,:),ic] = update_ic(gen,ic,best_dispatch,first_profile(t+1,:),dt(t),'initially constrained');%only updates the storage states in IC
    end
    for i = 1:1:n_g
        if ~(strcmp(gen(i).Type,'Chiller') || ismember(gen(i).Type,{'Electric Storage';'Thermal Storage';}))
            gen_output(2:end,i) = first_profile(2:end,i);
        end
    end
    dispatchable = logical(one_step.Organize.Dispatchable);
    gen_output = minimize_start_costs(gen,[date;forecast.Timestamp],dispatchable,gen_output,stor_power,alt_c,start_cost,dt,{'Chiller';},40);
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'Chiller') ||((ismember(gen(i).Type,{'Thermal Storage';}) && isfield(gen(i).QPform.output,'C')))
            first_profile(2:end,i) = gen_output(2:end,i);
        end
    end
    for i = 1:1:n_g
        if ismember(gen(i).Type,{'Thermal Storage';'Utility';'Electric Storage';'Solar';'Wind';})
            combinations_chill(:,i) = 1;  
        else
            combinations_chill(:,i) = first_profile(:,i)>0;
        end
    end
end
end %Ends function seperate_chiller_problem

function net_demand = agregate_demand(forecast)
net_demand = [];
if isfield(forecast,'Demand')
    outs = fieldnames(forecast.Demand);
    for j = 1:1:length(outs)
        net_demand.(outs{j}) = sum(forecast.Demand.(outs{j}),2);
    end
end
if isfield(forecast,'Building')
    outs2 = {'E';'H';'C';};
    nS = length(forecast.Timestamp);
    for j = 1:1:length(outs2)
        if ~isfield(net_demand,outs2{j})
            net_demand.(outs2{j}) = zeros(nS,1);
        end
        net_demand.(outs2{j}) = net_demand.(outs2{j}) + sum(forecast.Building.(strcat(outs2{j},'0')),2);
    end
end
end%Ends function agregate_demand

function v_h = value_heat(gen,dispatch,excess_heat,dt)
n_g = length(gen);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Thermal Storage') && strcmp(gen(i).Source,'Heat')
        loss = dt*(gen(i).QPform.Stor.SelfDischarge*gen(i).QPform.Stor.UsableSize);
        d_SOC = dispatch(2,i) - dispatch(1,i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
        if (d_SOC/dt + loss)<0 %discharging
            excess_heat = excess_heat -(d_SOC/dt + loss)*gen(i).QPform.Stor.DischEff;
        else %charging
            excess_heat = excess_heat -(d_SOC/dt +loss)/gen(i).QPform.Stor.ChargeEff; 
        end
    end
    if strcmp(gen(i).Type,'CHP Generator')
        j = 1;
        d = 0;
        states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
        if dispatch(2,i)>0
            excess_heat = excess_heat - gen(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
        end
        while j<=length(states) && dispatch(2,i)-d>0
            x_i = min(dispatch(2,i)-d,gen(i).QPform.(states{j}).ub(end));
            excess_heat = excess_heat + x_i*gen(i).QPform.output.H(min(j,length(gen(i).QPform.output.H(:,1))),1);
            d = d + x_i;
            j = j+1;
        end
    end
end
v_h = excess_heat<=1e-1;
end%Ends function value_heat

function stor_power = find_stor_power(gen,dispatch,dt)
n_g = length(gen);
stor_power = zeros(length(dt),n_g);  
for i = 1:1:n_g
    if isfield(gen(i).QPform,'Stor')
        for t = 1:1:length(dt)
            loss = dt(t)*(gen(i).QPform.Stor.SelfDischarge*gen(i).QPform.Stor.UsableSize);
            d_SOC = dispatch(t+1,i) - dispatch(t,i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
            if (d_SOC/dt(t) + loss)<0 %discharging
                stor_power(t,i) = -(d_SOC/dt(t) + loss)*gen(i).QPform.Stor.DischEff;
            else %charging
                stor_power(t,i) = -(d_SOC/dt(t) +loss)/gen(i).QPform.Stor.ChargeEff; 
            end
        end
    end
end
end%Ends function find_stor_power

function input = chill_input(qp_form,output)
net_out = 0;
input = 0;
i = 1;
while net_out<output
    seg = min(qp_form.(qp_form.states{i}).ub,output-net_out);
    net_out = net_out + seg;
    if isfield(qp_form.output,'H')
        input = input + seg*qp_form.output.H(i,2);
    end
    if isfield(qp_form.output,'E')
        if length(qp_form.output.E(1,:))>1
            input = input + seg*qp_form.output.E(i,2);
        else
            input = input + seg*qp_form.output.E;
        end
    end
    i = i+1;
end    
end%ends function gen_input