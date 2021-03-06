function qp = update_matrices_step(gen,building,cool_tower,subnet,options,qp,forecast,scale_cost,marginal,stor_power,dt,ic,first_profile,limit,t,temperatures)
% update the equalities with the correct demand, and scale fuel and electric costs
% ec is the expected end condition at this time stamp (can be empty)
% stor_power is the expected output/input of any energy storage at this timestep (can be empty)
% min_power and MaxPower define the range of this generator at this timestep
% temperatures is the building and cooling tower water loop temperatures
qp.solver = options.solver;
n_g = length(gen);
n_b = length(building);
n_l = length(qp.Organize.Transmission);
n_ct = length(cool_tower);

upper_bound = zeros(1,n_g);
dx_dt = zeros(1,n_g);
for i = 1:1:n_g
    if ~isempty(gen(i).QPform.states)
        if strcmp(gen(i).Type,'Hydro Storage')
            %%do something to set UB
        else
            states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,end))),end);
            for j = 1:1:length(states)
                upper_bound(i) = upper_bound(i) + gen(i).QPform.(states{j}).ub(end);
            end
        end
    end
    if isfield(gen(i).VariableStruct,'dX_dt')
        dx_dt(i) = gen(i).VariableStruct.dX_dt;
    end
end
if strcmp(limit, 'unconstrained')
    min_power = 0*upper_bound;
    max_power = upper_bound;
elseif strcmp(limit, 'constrained')
    min_power = max(0,ic(1:n_g)-dx_dt*dt(t));
    max_power = min(upper_bound,ic(1:n_g)+dx_dt*dt(t));
else
    min_power = max(0,ic(1:n_g)-dx_dt*sum(dt(1:t)));
    max_power = min(upper_bound,ic(1:n_g)+dx_dt*sum(dt(1:t)));
end

network_names = fieldnames(subnet);
%% update demands and storage self-discharge
for net = 1:1:length(network_names)
    out = subnet.(network_names{net}).abbreviation;
    if strcmp(network_names{net},'Hydro')
        %don't do a water balance, since it depends on multiple time steps.
        %Any extra outflow at this time step is subtracted from the expected
        %outflow at the next step (same SOC and flows up river).
    else
        if strcmp(network_names{net},'Electrical') && options.SpinReserve
            qp.b(qp.Organize.SpinReserve) = -options.SpinReservePerc/100*sum(forecast.Demand.E(t,:),2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        end
        qp.constDemand.(out).req = zeros(length(subnet.(network_names{net}).nodes),1);
        qp.constDemand.(out).load = zeros(length(subnet.(network_names{net}).nodes),n_g);
        for n = 1:1:length(subnet.(network_names{net}).nodes)
            equip = subnet.(network_names{net}).Equipment{n};
            req = qp.Organize.Balance.(network_names{net})(n);
            loads = nonzeros(qp.Organize.Demand.(network_names{net}){n});
            if ~isempty(loads)%need this to avoid the next line when there is no load
                qp.beq(req) = sum(forecast.Demand.(out)(t,loads),2); %multiple demands can be at the same node
            end
            for j = 1:1:length(equip)
                k = equip(j);
                if strcmp(network_names{net},'Electrical') 
                    if strcmp(gen(k).Source,'Renewable')% subtract renewable generation 
                        qp.beq(req) = qp.beq(req) - forecast.Renewable(t,k); %put renewable generation into energy balance at correct node
                    end
                end
                if ~isempty(first_profile) && ismember(gen(k).Type,{'Electric Storage';'Thermal Storage';})
                    if isfield(gen(k).QPform.output,out)
                        loss = (gen(k).QPform.Stor.SelfDischarge*gen(k).QPform.Stor.UsableSize*gen(k).QPform.Stor.DischEff);
                        d_soc = first_profile(t+1,k) - ic(k);%positive means charging
                        qp.beq(req) = qp.beq(req) + (d_soc/dt(t)+loss)*gen(k).QPform.Stor.DischEff;
                        qp.b(qp.Organize.Inequalities(k)) = -(1/gen(k).QPform.Stor.ChargeEff-gen(k).QPform.Stor.DischEff)*(d_soc/dt(t)+loss);
                    end
                end
                if isfield(gen(k).QPform,'constDemand') && isfield(gen(k).QPform.constDemand,out)
                    qp.constDemand.(out).req(n,1) = req;
                    qp.constDemand.(out).load(n,k) = gen(k).QPform.constDemand.(out);
                end
            end
        end
    end
end  

%% Building inequalities & equality
for i = 1:1:n_b
    states = qp.Organize.States{n_g+n_l+i};
    req = qp.Organize.Building.Electrical.req(i);
    qp.beq(req) = qp.beq(req) + forecast.Building.E0(t,i);%Equipment and nominal Fan Power
    %Heating Equality
    req = qp.Organize.Building.DistrictHeat.req(i);
    qp.Organize.Building.H_Offset(1,i) = building(i).QPform.UA*(forecast.Building.Tzone(t,i)-forecast.Building.Tset_H(t,i)) - forecast.Building.H0(t,i);
    qp.beq(req) = qp.beq(req) - qp.Organize.Building.H_Offset(1,i);%Heating load
    qp.lb(states(2)) = qp.lb(states(2)) + qp.Organize.Building.H_Offset(1,i);
    
    %Cooling Equality
    req = qp.Organize.Building.DistrictCool.req(i);
    qp.Organize.Building.C_Offset(1,i) = building(i).QPform.UA*(forecast.Building.Tset_C(t,i)-forecast.Building.Tzone(t,i)) - forecast.Building.C0(t,i);
    qp.beq(req) = qp.beq(req) - qp.Organize.Building.C_Offset(1,i);%Cooling Load
    qp.lb(states(3)) = qp.lb(states(3)) + qp.Organize.Building.C_Offset(1,i);
    
    %heating inequality
    qp.b(qp.Organize.Building.r(i)) = building(i).QPform.UA*forecast.Building.Tset_H(t,i) + building(i).QPform.Cap*temperatures.build(1,i)/(3600*dt(t));
    qp.A(qp.Organize.Building.r(i),states(1)) = (building(i).QPform.UA+building(i).QPform.Cap/(3600*dt(t)));
    %cooling inequality
    qp.b(qp.Organize.Building.r(i)+1) = -building(i).QPform.UA*forecast.Building.Tset_C(t,i) + building(i).QPform.Cap*temperatures.build(1,i)/(3600*dt(t));
    qp.A(qp.Organize.Building.r(i)+1,states(1)) = -building(i).QPform.UA - building(i).QPform.Cap/(3600*dt(t));
    %penalty states (excess and under temperature)
    qp.b(qp.Organize.Building.r(i)+2) = forecast.Building.Tmax(t,i);%upper buffer inequality, the longer the time step the larger the penaly cost on the state is.
    qp.b(qp.Organize.Building.r(i)+3) = -forecast.Building.Tmin(t,i);%lower buffer inequality
end
%% Cooling Tower water loop equality
for i = 1:1:n_ct
    %1st order temperature model: 0 = -T(k) + T(k-1) + dt/Capacitance*(energy balance)
    state = qp.Organize.States{n_g+n_l+n_b+i};
    req = qp.Organize.Balance.CoolingWater(i);
    qp.Organize.cool_tower.req(i) = req;
    capacitance = cool_tower(i).fluid_capacity*cool_tower(i).fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
    qp.Aeq(req,:) = qp.Aeq(req,:)*dt(t)/capacitance; %energy balance for chillers & cooling tower fans already put in this row of Aeq
    qp.Aeq(req,state) = -1; %-T(k)
    qp.beq(req) = -temperatures.cool_tower(i);%-T(k-1)
end

%% Update upper and lower bounds based on ramping constraint (if applicable)
if ~isempty(first_profile)
    for i = 1:1:n_g
        states = qp.Organize.States{i};
        if ismember(gen(i).Type,{'Utility';'Electric Generator';'CHP Generator';'Chiller';'Heater';})%all generators and utilities
            if min_power(i)>sum(qp.lb(states))
                qp.Organize.Dispatchable(i) = 0; %can't shut off
                %raise the lower bounds
                power_add = min_power(i) - sum(qp.lb(states));
                for f = 1:1:length(states) %starting with lb of 1st state in generator, then moving on
                    gap = qp.ub(states(f)) - qp.lb(states(f));
                    if gap>0
                        gap = min(gap,power_add);
                        qp.lb(states(f)) = qp.lb(states(f)) + gap;
                        power_add = power_add-gap;
                    end
                end
            end
            if max_power(i)<gen(i).Size %lower the upper bounds
                power_sub = sum(qp.ub(states)) - max_power(i);
                for f = length(states):-1:1 %starting with lb of 1st state in generator, then moving on
                    if qp.ub(states(f))>0
                        gap = min(qp.ub(states(f)),power_sub);
                        qp.ub(states(f)) = qp.ub(states(f)) - gap;
                        power_sub = power_sub-gap;
                    end
                end
            end
        elseif any(strcmp(gen(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';}))
            %update storage output range to account for what is already scheduled
            if isfield(gen(i).QPform.output,'H')
                req = qp.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            elseif isfield(gen(i).QPform.output,'E')
                req = qp.Organize.Balance.Electrical;
            elseif isfield(gen(i).QPform.output,'C')
                req = qp.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
            end
            s = states(1);
            qp.lb(s) = qp.lb(s) - stor_power(i);%change in storage for this power output
            qp.ub(s) = qp.ub(s) - stor_power(i);%change in storage for this power output
            chargingSpace = max(0,sum((gen(i).QPform.Stor.UsableSize-first_profile(t+1,i)).*qp.Aeq(req,s))/dt(t));%you can't charge more than you have space for
            qp.lb(s) = max(qp.lb(s), -chargingSpace);
            qp.ub(s) = min(qp.ub(s), first_profile(t+1,i)/dt(t));%you can't discharge more than you have stored            
            if qp.Organize.SpinRow(i)~=0 %update spinning reserve (max additional power from storage
                qp.beq(qp.Organize.SpinRow{i}) = qp.ub(s);
            end
        end
        if strcmp(gen(i).Type,{'Hydro Storage';})
            %update the equality with NewPower*conversion + spill - Outflow = nominal PowerGen Flow
            qp.beq(qp.Organize.HydroEqualities(i)) = -first_profile(i)*gen(i).QPform.Stor.Power2Flow;
        end
        qp.lb(states) = min(qp.lb(states),qp.ub(states)); %fix rounding errors
    end
end

%% Adding cost for Line Transfer Penalties to differentiate from spillway flow 
%%also want to penalize spill flow to conserve water for later, spill flow
%%is just to meet instream requirments when the power gen needs to be low
if ismember('Hydro',network_names) && ~isempty(subnet.Electrical.lineNames)
    for k = 1:1:length(subnet.Electrical.lineNames)
        line = subnet.Electrical.lineNumber(k);
        states = qp.Organize.States{n_g+line};
        if length(states)>1 %bi-directional transfer with penalty states
            s2 = states(2):qp.Organize.t1States:(states(2));%penalty from a to b
            s3 = states(3):qp.Organize.t1States:(states(3));% penalty from b to a
            qp.f(s2) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
            qp.f(s3) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
        end
    end 
    for i = 1:1:n_g
        if strcmp(gen(i).Type,{'Hydro Storage';}) && isfield(gen(i).QPform,'S') 
            states = qp.Organize.States{i};
            qp.f(states(2)) = 0.0001/gen(i).QPform.Stor.Power2Flow;
        end
    end
end 

%% update costs
H = diag(qp.H);%convert to vector
for i = 1:1:n_g
    if ~isempty(qp.Organize.States{i})
        states = qp.Organize.States{i}; %states of this generator
        if strcmp(gen(i).Type,'Utility') && strcmp(gen(i).Source,'Electricity') 
            qp.f(states(1)) = qp.f(states(1))*scale_cost(i)*dt(t);
            if gen(i).VariableStruct.SellBackRate == -1
                qp.f(states(2)) = qp.f(states(2))*scale_cost(i)*dt(t); %sellback is a fixed percent of purchase costs (percent set wehn building matrices)
            else
                %constant sellback rate was taken care of when building the matrices
            end
        elseif ismember(gen(i).Type,{'Electric Generator';'CHP Generator';'Utility';'Chiller';'Heater';})%all generators and utilities
            H(states) = H(states)*scale_cost(i)*dt(t);
            qp.f(states) = qp.f(states)*scale_cost(i)*dt(t);
        elseif ismember(gen(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})%all  storage
            stor_cat = fieldnames(gen(i).QPform.output);
            %penalize deviations from the expected storage output
            %don't penalize spinning reserve state
            if isempty(first_profile)
                qp.f(states(1)) = marginal.(stor_cat{1});
                if any(strcmp(network_names,'Electrical')) && strcmp(stor_cat{1},'H') % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                    H(states(1)) = 0;
                elseif strcmp(stor_cat{1},'W') %Hydro needs to have high cost, but not too high that it never uses it; set arbitrary number
                    H(states(1)) = 1e5*marginal.(stor_cat{1});
                else 
                    H(states(1)) = 1e8*marginal.(stor_cat{1});
                end
            else
                a = 4; %sets severity of quadratic penalty
                MaxCharge = gen(i).QPform.Ramp.b(1); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                if gen(i).QPform.Stor.UsableSize>first_profile(t+1,i)
                    MaxCharge = min((gen(i).QPform.Stor.UsableSize-first_profile(t+1,i))/dt(t),MaxCharge);
                end
                qp.f(states(1)) = marginal.(stor_cat{1});
                if MaxCharge == 0
                    H(states(1)) = 0;
                else
                    H(states(1)) = a*2*qp.f(states(1))/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
                qp.f(states(2)) = 1e-6; %small penalty on charging state so that if there are no other generators it keeps the charging constraint as an equality
            end
        end
    end
end 
qp.constCost = qp.constCost.*scale_cost*dt(t);
if options.SpinReserve
    sr_short = qp.Organize.SpinReserveStates(1,n_g+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    sr_ancillary = qp.Organize.SpinReserveStates(1,n_g+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if options.SpinReservePerc>5 %more than 5% spinning reserve
        spin_cost = 2*dt(t)./(options.SpinReservePerc/100*sum(forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        spin_cost = 2*dt(t)./(0.05*sum(forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(sr_short) = spin_cost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    qp.f(sr_short) = 0.05*dt(t); % $0.05 per kWh
end
qp.H = diag(H);%convert back to matrix
end%ends function update_matrices_step