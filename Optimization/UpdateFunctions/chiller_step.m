function qp = chiller_step(gen,building,cool_tower,subnet,sr,qp,marginal,first_profile)
% update the equalities with the correct demand, and scale fuel and electric costs
% ec is the expected end condition at this time stamp (can be empty)
% stor_power is the expected output/input of any energy storage at this timestep (can be empty)
% min_power and MaxPower define the range of this generator at this timestep
n_g = length(gen);
ab_chiller = [];
sub_problem = true(1,n_g);
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Chiller')
        states = qp.Organize.States{i};
        if isfield(gen(i).QPform.output,'H')
            qp.f(states) = -marginal.H*gen(i).QPform.output.H(:,2);
            ab_chiller(end+1) = i;
        else
            qp.f(states) = -marginal.E*gen(i).QPform.output.E(:,2);
        end
    elseif strcmp(gen(i).Type,'Cooling Tower')
        states = qp.Organize.States{i};
        qp.f(states) = -marginal.E*gen(i).QPform.output.E(:,2);
    elseif strcmp(gen(i).Type,'Thermal Storage') && isfield(gen(i).QPform.output,'C')
        states = qp.Organize.States{i};
        %included but no costs
    else
        sub_problem(i) = false;
    end
    %%Do I need to add something to keep the cooling part of a building, or
    %%is that kept automatically?
end

%% set upper limit on ab chiller
if ~isempty(ab_chiller)
    req = qp.Organize.Balance.DistrictHeat;
    excess_heat = -qp.beq(req);
    for i = 1:1:n_g
        if ismember(gen(i).Type,{'CHP Generator';'Heater'}) && first_profile(i)>0
            states = qp.Organize.States{i};
            D = 0;
            j = 1;
            if isfield(gen(i).QPform,'constDemand')
                excess_heat = excess_heat - gen(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
            end
            while D<first_profile(i) && j<length(states)
                g = min(first_profile(i)-D,qp.ub(states(j)));
                D = D+g;
                excess_heat = excess_heat + g*qp.Aeq(req,states(j));
                j = j+1;
            end
        elseif strcmp(gen(i).Type,'Chiller')
            if isfield(gen(i).QPform.constDemand,'H')
                qp.constCost(i) = qp.constCost(i) + gen(i).QPform.constDemand.H*marginal.H;
            end
            if isfield(gen(i).QPform.constDemand,'E')
                qp.constCost(i) = qp.constCost(i) + gen(i).QPform.constDemand.E*marginal.E;
            end
        end
    end
    %solve for max abChill output for excess heat
    for k = 1:1:length(ab_chiller)
        i = ab_chiller(k);
        states = qp.Organize.States{i};
        for j=1:1:length(states)
            if j == 1
                if excess_heat>(gen(i).QPform.constDemand.H+qp.lb(states(1))*(-gen(i).QPform.output.H(1,2)))%enough heat to get above LB
                    excess_heat = excess_heat - gen(i).QPform.constDemand.H;
                    qp.constCost(i) = 0;
                else
                    excess_heat = 0;
                end
            end
            if qp.ub(states(j))>0
                H = min(excess_heat,qp.ub(states(j))*(-gen(i).QPform.output.H(j,2)));
                excess_heat = excess_heat - H;
                free_heat = H/(-gen(i).QPform.output.H(j,2));
                qp.f(states(j)) = (1-free_heat/qp.ub(states(j)))*qp.f(states(j));%set cost to zero if there is 'free excess heat', but don't reduce the available capacity
            end
        end
    end
end
qp.Organize.HeatVented = [];
%remove all non chillers and cold storage & add cost
qp = rmfield(qp,'constDemand');
qp = disable_generators(qp,[],sub_problem);

%eliminate non district cooling energy balances & other constraints
r = length(qp.b);
x_keep = true(length(qp.f),1);
req_keep = false(length(qp.beq),1);
r_keep = true(r,1);
req_keep(qp.Organize.Balance.DistrictCool) = true;
if qp.excessHeat
    x_keep(qp.Organize.HeatVented) = false;
end
%%remove transmission states and constraints
n_b = length(building);
n_l = length(qp.Organize.Transmission);
for i = 1:1:n_l
    if ~ismember(i,subnet.DistrictCool.lineNumber)
        x_keep(qp.Organize.States{n_g+i}) = false;
        r_keep(qp.Organize.Transmission(i)) = false;
        r_keep(qp.Organize.Transmission(i)+1) = false;
    end
end
if isfield(subnet,'Electrical') && sr
    r_keep(qp.Organize.SpinReserve) = false;
end
for i = 1:1:n_b
    r_keep(qp.Organize.Building.r(i)) = false;
end   
qp.H = qp.H(x_keep,x_keep);
qp.f = qp.f(x_keep);
qp.Aeq = qp.Aeq(req_keep,x_keep);
qp.beq = qp.beq(req_keep);
if r>0 && nnz(r_keep)>0
    qp.A = qp.A(r_keep,x_keep);
    qp.b = qp.b(r_keep);
else qp.A = [];
    qp.b = [];
end
qp.lb = qp.lb(x_keep);
qp.ub = qp.ub(x_keep);
qp.Organize.t1States = length(qp.f);
qp.Organize.t1Balances = length(subnet.DistrictCool.nodes);
qp.Organize.Balance = [];
qp.Organize.Balance.DistrictCool = 1;
