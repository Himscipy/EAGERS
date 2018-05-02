function combinations = create_combinations(gen,options,qp,net_demand,ic,dt,t,k_c,out)
%% test if each combination is capable of meeting demand
%if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
%if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
%%could improve this section to check if feasible with network losses (not sure how)
if strcmp(out,'C')
    c = net_demand.C;
    net_demand = [];
    net_demand.C = c;
    include = {'Chiller'};
    req = qp.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    [combinations,~] = build_cases(gen,include,qp.Organize.Dispatchable,qp.Organize.Enabled,[]);
    [combinations,~] = check_feas(gen,options,combinations,qp,req,net_demand.C(t),[],ic,dt);
elseif strcmp(out,'E')
    if isfield(net_demand,'E') 
        net_demand.E = net_demand.E;
    else 
        net_demand.E = 0;
    end
    if isfield(net_demand,'DC')
        net_demand.E = net_demand.E + net_demand.DC;
    end
    include = {'Electric Generator';'CHP Generator';'Hydrogen Generator';'Heater';'Electrolyzer'};%no utilities or storage systems are dispatchable
    if isempty(k_c)
        include(end+1:end+3) = {'Chiller';'Absorption Chiller';'Cooling Tower';};
    else
        net_demand = rmfield(net_demand,'C');
    end
    req = qp.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    
    n_g = length(gen);
    is_chp = false;
    for i = 1:1:n_g
        if strcmp(gen(i).Type,'CHP Generator')
            is_chp = true;
        end
    end
    [all_combinations,~] = build_cases(gen,include,qp.Organize.Dispatchable,qp.Organize.Enabled,[]);
    [all_combinations,~] = check_feas(gen,options,all_combinations,qp,req,net_demand.E(t),[],ic,dt);
    if isempty(k_c)
       combinations = all_combinations;
       if is_chp && ~isempty(combinations)
            req = qp.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            [combinations,~] = check_feas(gen,options,combinations,qp,req,net_demand.H(t),'H',ic,dt);
       end
       if isfield(net_demand,'C') && ~isempty(combinations)
           req = qp.Organize.Balance.DistrictCool; %rows of Aeq associated with cooling demand
           [combinations,~] = check_feas(gen,options,combinations,qp,req,net_demand.C(t),'C',ic,dt);
       end
    else
        [combinations,inc] = build_cases(gen,include,qp.Organize.Dispatchable,qp.Organize.Enabled,k_c);
        [combinations,notFeas] = check_feas(gen,options,combinations,qp,req,net_demand.E(t),[],ic,dt);
        if isempty(combinations) && ~isempty(notFeas)
            for i = 1:1:length(notFeas(:,1))%%need extra logic in this create combinations so that a set of power/heat generators that is not feasible with the ideal chiller set may be feasible with a different chiller set that also meets the cooling constraints
                if any(notFeas(i,inc'))
                    f_rows = ismember(all_combinations(:,inc'),notFeas(i,inc'),'rows');%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                else
                    f_rows = all(~all_combinations(:,inc'),2);%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                end
                combinations(end+1:end+nnz(f_rows),:) = all_combinations(f_rows,:);%add these rows to K
            end
        end
        if is_chp && ~isempty(combinations)
            req = qp.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            [combinations,notFeas] = check_feas(gen,options,combinations,qp,req,net_demand.H(t),'H',ic,dt);
            if isempty(combinations) && ~isempty(notFeas) && ~isempty(all_combinations)
                for i = 1:1:length(notFeas(:,1))
                    if any(notFeas(i,inc'))
                        f_rows = ismember(all_combinations(:,inc'),notFeas(i,inc'),'rows');%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                    else
                        f_rows = all(~all_combinations(:,inc'),2);%find any rows of K_all with the infeasible combo of electric & CHP & heaters
                    end
                    combinations(end+1:end+nnz(f_rows),:) = all_combinations(f_rows,:);%add these rows to K
                end
            end
        end
    end
end
if isempty(combinations)
    disp('Zero feasible combinations to test: ERROR')
end
%All combinations at this point should be feasible

% Premise (Maybe?), any extra generators beyond what is necessary to meet demand, or
% that are more expensive than grid are unneccessary
%if a combination is feasible, and another combination includes all those
%generators + additional  generators, those other cases can be removed
end%Ends function create_combinations

function [combinations,not_feas] = check_feas(gen,options,combinations,qp,req,net_demand,out,ic,dt)
n_g = length(combinations(1,:));
limit_lower = zeros(1,n_g);
limit_upper = zeros(1,length(limit_lower));
for i = 1:1:length(limit_lower)
    states = qp.Organize.States{i};
    if ~isempty(states) && any(any(qp.Aeq(req,states)~=0))
        row_i = find(any(qp.Aeq(req,states)~=0));
        if strcmp(gen(i).Type,'Utility')%if it is a utility
            if length(states) ==2
                limit_lower(i) = sum(qp.Aeq(req(row_i),states(2))*qp.ub(states(2)));%if there is sellback, take 2nd state
                limit_upper(i) = sum(qp.Aeq(req(row_i),states(1))*qp.ub(states(1)));
            else
                limit_lower(i) = sum(qp.Aeq(req(row_i),states).*qp.lb(states)');
                limit_upper(i) = sum(qp.Aeq(req(row_i),states).*qp.ub(states)');
            end
        elseif any(strcmp(gen(i).Type,{'Hydro Storage';'Electric Storage';'Thermal Storage';}))%if it is storage
            s = states(1);
            limit_lower(i) = sum(qp.Aeq(req(row_i),s)*qp.lb(s));
            limit_upper(i) = sum(qp.Aeq(req(row_i),s)*qp.ub(s));
            if ~isempty(ic)%bounds reduced by state of charge
                limit_upper(i) = min(ic(i).*qp.Aeq(req,s)/dt,limit_upper(i));
                chargingSpace = sum((gen(i).QPform.Stor.UsableSize-ic(i)).*qp.Aeq(req,s));
                limit_lower(i) = max(-chargingSpace/dt, limit_lower(i));
            end
        else
            limit_lower(i) = sum(qp.Aeq(req(row_i),states).*qp.lb(states)');
            limit_upper(i) = sum(qp.Aeq(req(row_i),states).*qp.ub(states)');
            if limit_upper(i)<limit_lower(i)%if it consumes this demand (i.e. a chiller when examining electrical demand)
                temp = limit_upper(i);
                limit_upper(i) = limit_lower(i);
                limit_lower(i) = temp;
            end
        end
        if isfield(gen(i).QPform,'constDemand') && ~isempty(out) && isfield(gen(i).QPform.constDemand,out)
            limit_upper(i) = limit_upper(i) - gen(i).QPform.constDemand.(out);
            limit_lower(i) = limit_lower(i) - gen(i).QPform.constDemand.(out);
        end
    end
end
if ~isempty(out) && strcmp(out,'H') && options.excessHeat
    limit_lower = limit_lower - inf;
end
n = length(combinations(:,1));
a = (ones(n,1)*limit_upper);
a(combinations==0) = 0;
sum_upper_bound = sum(a,2);
b = (ones(n,1)*limit_lower);
b(combinations==0) = 0;
sum_lower_bound = sum(b,2);

% %% Trying to see if we can add transmission losses
% if ~isempty(QP.Organize.Transmission)
%     [sumUB,sumLB] = TransferLoss(QP,A,B,sumUB,sumLB);
% end 

keep = ((sum_upper_bound>=net_demand) & (sum_lower_bound<=net_demand));%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
not_feas = combinations(~keep,:);
combinations = combinations(keep,:);
end%ends function check_feas

function [combinations,inc] = build_cases(gen,include,dispatchable,available,k_c)
n_g = length(gen);
inc = [];
ninc = [];
for i = 1:1:n_g
    if ismember(gen(i).Type,include) && dispatchable(i) && available(i)
        inc(end+1) = i;
    elseif available(i)
        ninc(end+1) = i;
    end
end
n = length(inc);
combinations = zeros(2^n,n_g); %K is a matrix of all the possible generator on/off combinations 
if ~isempty(ninc)
    for j = 1:1:length(ninc)
        if ~isempty(k_c) && strcmp(gen(ninc(j)).Type,'Chiller')
            if k_c(ninc(j))
                combinations(:,ninc(j)) = ninc(j); % Chiller status determined earlier
            end
        else
            combinations(:,ninc(j)) = ninc(j); % all systems that are always included
        end
    end
end

for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=2^n
        combinations(r+1:r+z,inc(j)) = inc(j);
        r=r+2*z;
    end
end
end%ends function build_cases