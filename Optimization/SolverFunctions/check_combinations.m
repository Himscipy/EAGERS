function [best_dispatch,alt] = check_combinations(qp_0,combinations,alt,ic,net_demand,dt,t)
%This function identifies all feasible generator combinations that can meet
%the demand at this moment in time
%These feasible combinations are then tested, begining with the options
%that have the fewest number of generators
%Some combinations are avoided if they can be pre-emptively determined to
%be more costly
%The function returns the feasible generator dispatches, the cost of that 
%dispatch, and the on/off binary matrix that represents those combinations
[lines,~] = size(combinations);
if lines>1 && license('test','Distrib_Computing_Toolbox') 
    parallel = true;
else
    parallel = false;
end

n_g = length(qp_0.constCost);
[~,n] = size(qp_0.organize);
n_b = length(qp_0.Organize.Building.r);
n_l = length(qp_0.Organize.Transmission);
n_h = nnz(qp_0.Organize.Hydro);

dispatch = zeros(lines,n);
line_loss = zeros(lines,n_l);
excess_heat = zeros(lines,nnz(qp_0.Organize.HeatVented));
excess_cool = zeros(lines,nnz(qp_0.Organize.CoolVented));
hydro_soc = zeros(lines,n_h);
temperature = zeros(lines,n_b);
heating = zeros(lines,n_b);
cooling = zeros(lines,n_b);
cost = zeros(lines,1);
feasible = false(lines,1);
flag1 = zeros(lines,1);

if parallel
    parfor i = 1:lines
        QP = disable_generators(qp_0,[],combinations(i,:));%Disable generators here
        [x, flag1(i)] = call_solver(QP);
        if flag1(i)==1
            cost(i) = 0.5*x'*QP.H*x + x'*QP.f;
            [dispatch(i,:),line_loss(i,:),excess_heat(i,:),excess_cool(i,:),hydro_soc(i,:),temperature(i,:),heating(i,:),cooling(i,:)] = sort_solution_step(x,QP);
        end
    end
    cost = cost+sum(ones(lines,1)*qp_0.constCost.*(combinations>0),2);
else
    nzK = sum(combinations>0,2);%this is the number of active generators per combination (nonzeros of K)
    [~, line] = sort(nzK); %sort the rows by number of generators that are on
    combinations = combinations(line,:);
    i = 1;
    while i<=length(line)
        QP = disable_generators(qp_0,[],combinations(i,:));%Disable generators here
        [x, flag1(i)] = call_solver(QP);
        if flag1(i)==1
            cost(i) = 0.5*x'*QP.H*x + x'*QP.f;
            cost(i) = cost(i)+sum(qp_0.constCost.*(combinations(i,:)>0),2);
            [dispatch(i,:),line_loss(i,:),excess_heat(i,:),excess_cool(i,:),hydro_soc(i,:),temperature(i,:),heating(i,:),cooling(i,:)] = sort_solution_step(x,QP);
        end
        if i<length(line)
            combinations = eliminate_cases(combinations,QP,net_demand,i,cost(i),min(cost(1:i-1)),dt,t);
        end
        i = i+1;
    end
end

feasible(flag1==1) = true;
cost(feasible==false) = inf;
[cost,I] = sort(cost);
nP = nnz(feasible==true);
combinations = combinations(I(1:nP),:); %the n best combinations

alt.Binary{t} = combinations>0;
alt.Disp{t} = dispatch(I(1:nP),1:n_g);
alt.Cost{t} = cost(1:nP)-cost(1);
alt.Generators{t} = dispatch(I(1:nP),1:n_g);
alt.LineFlows{t} = dispatch(I(1:nP),n_g+1:n_g+n_l);
alt.Buildings{t} = dispatch(I(1:nP),n_g+n_l+1:n_g+n_l+n_b);
alt.LineLoss{t} = line_loss(I(1:nP),:);
alt.excessHeat{t} = excess_heat(I(1:nP),:);
alt.excessCool{t} = excess_cool(I(1:nP),:);
alt.hydroSOC{t} = hydro_soc(I(1:nP),:);
alt.Temperature{t} = temperature(I(1:nP),:);
alt.Heating{t} = heating(I(1:nP),:);
alt.Cooling{t} = cooling(I(1:nP),:);
if isempty(alt.Disp{t})
    disp(['No feasible combination of generators at step' num2str(t)]);
    best_dispatch = ic;
else
    best_dispatch = alt.Disp{t}(1,:);
end
end %ends function check_combinations

function combinations = eliminate_cases(combinations,QP,net_demand,i,cost,best_cost,dt,t)
%%reduce size of K if you don't have parallel computing toolbox
if isempty(best_cost)
    best_cost = inf;
end
if cost < best_cost
    %remove combinations that are the best combination and additional more expensive generators
    outs = fieldnames(net_demand);
    n_g = length(QP.constCost);
    for s = 1:1:length(outs)
        req = [];
        if strcmp(outs{s},'E')
            req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
        elseif strcmp(outs{s},'H')
            req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
        elseif strcmp(outs{s},'C')
            req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
        elseif strcmp(outs{s},'W')
            req = QP.Organize.Balance.Hydro; %rows of Aeq associated with hydro demand
        end

        if length(req) ==1 % The following definitely works with only 1 node
            best_cost_per_kWh = cost/(net_demand.(outs{s})(t)*dt);
            n_row = length(combinations(:,1))-i;%how many rows do you have left
            include = false(1,n_g);
            min_rate = inf+zeros(1,n_g);
            for j = 1:1:n_g
                if QP.Organize.Dispatchable(j)
                    states = QP.Organize.States{j};
                    if any(QP.Aeq(req,states))
                        include(j) = true;
                        min_rate(j) = QP.f(states(1));
                        if min_rate(j)== 0 && QP.Aeq(QP.Organize.Balance.Electrical,states(1))<0
                            min_rate(j) = -.1*QP.Aeq(QP.Organize.Balance.Electrical,states(1));
                        elseif min_rate(j)== 0 && QP.Aeq(QP.Organize.Balance.DistrictHeat,states(1))<0
                            min_rate(j) = -.1*QP.Aeq(QP.Organize.Balance.DistrictHeat,states(1));
                        end
                    end
                end
            end
            if i<length(combinations(:,1))
                pos_cheaper_gen = nonzeros((1:n_g).*include.*((min_rate<best_cost_per_kWh) & (combinations(i,:)==0)));%possibly cheaper generators that are not on for this case, but could be on
            else
                pos_cheaper_gen = [];
            end
            if ~isempty(pos_cheaper_gen)
                ck_rows = ~any((ones(n_row,1)*combinations(i,include)-combinations(i+1:end,include))>0,2); %identify future rows that have the current set of active generators + more
                no_cheaper = ~any(combinations(i+1:end,pos_cheaper_gen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
                rm_rows = ck_rows & no_cheaper;
                if nnz(rm_rows)>0
                    combinations = combinations([true(i,1);~rm_rows],:);  
                end
            end
            %remove combinations that swap the current
            %combination with a more expensive generator
            for m = 1:1:n_g
                if include(m) && combinations(i,m)>0 && i<length(combinations(:,1))
                    expensive_gens = nonzeros((1:n_g).*include.*(min_rate>min_rate(m)));
                    for q = 1:1:length(expensive_gens)
                        if combinations(i,expensive_gens(q))==0
                            add_gen = zeros(1,n_g);
                            add_gen(expensive_gens(q)) = expensive_gens(q);
                            add_gen(m) = -m;
                            swap_combo = combinations(i,:)+add_gen;
                            keep = ~ismember(combinations,swap_combo,'rows');
                            combinations = combinations(keep,:);
                        end
                    end
                end
            end
        end
    end
end
end %end function eliminate_cases