function solution = sort_solution(x,qp)
[m,n] = size(qp.organize);
n_s = m-1;
n_g = length(qp.constCost(1,:));
n_b = length(qp.Organize.Building.r);
n_l = length(qp.Organize.Transmission);
n_h = nnz(qp.Organize.Hydro);
solution.Dispatch = zeros(n_s+1,n_g);
solution.hydroSOC = zeros(n_s,n_h);
solution.excessHeat = zeros(n_s,nnz(qp.Organize.HeatVented));
solution.excessCool = zeros(n_s,nnz(qp.Organize.CoolVented));
solution.LineFlows = zeros(n_s,n_l);
solution.LineLoss = zeros(n_s,n_l);
solution.Buildings.Heating = zeros(n_s,n_b);
solution.Buildings.Cooling = zeros(n_s,n_b);
solution.Buildings.Temperature = zeros(n_s,n_b);
for i = 1:1:n_g
    if isfield(qp,'Renewable') && any(qp.Renewable(:,i)~=0)
        solution.Dispatch(2:end,i) = qp.Renewable(:,i);
    else
        out_vs_state = qp.Organize.Out_vs_State{1,i};%linear maping between state and output
        for t = 1:1:n_s+1
            if ~isempty(qp.organize{t,i})
                states = qp.organize{t,i};
                p = 0;
                for j = 1:1:length(states)
                    p = p + out_vs_state(j)*x(states(j)); %record this combination of outputs (SOC for storage)
                end
                if abs(p)>2e-4%added this to avoid rounding errors in optimization. When generator is locked off, the UB is set to 1e-5
                    solution.Dispatch(t,i) = p;
                end
            end
        end
    end
end
for i = 1:1:length(qp.Organize.Hydro)%Get SOC of each generator into a matrix for all time steps
    for t = 1:1:n_s
        solution.hydroSOC(t,i) = x(qp.organize{t+1,qp.Organize.Hydro(i)}+1,1)+qp.Organize.hydroSOCoffset(i);
    end
end
for i = 1:1:n_l
    for t = 1:1:n_s
        solution.LineFlows(t,i) = sum(x(qp.organize{t+1,i+n_g}));%power transfer
        if length(qp.Organize.States{i+n_g})>1
            solution.LineLoss(t,i) = sum(x(qp.organize{t+1,i+n_g}+1)); %down (positive) lines
            solution.LineLoss(t,i) = solution.LineLoss(t,i) + sum(x(qp.organize{t+1,i+n_g}+2)); %up (negative) lines
        end
    end
end
for i = 1:1:n_b
    states = qp.Organize.States{n_g+n_l+i};
    temperature_set = states(1):qp.Organize.t1States:(n_s-1)*qp.Organize.t1States + states(1);
    for t = 1:1:n_s
        solution.Buildings.Temperature(:,i) = x(temperature_set,1);
        solution.Buildings.Heating(:,i) = x(temperature_set+1,1) - qp.Organize.Building.H_Offset(:,i);
        solution.Buildings.Cooling(:,i) = x(temperature_set+2,1) - qp.Organize.Building.C_Offset(:,i);
    end
end

solution.Dispatch(abs(solution.Dispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
%pull out any dumped heat
for i = 1:1:length(qp.Organize.HeatVented)
    if qp.Organize.HeatVented(i)>0
        for t = 1:1:n_s
            solution.excessHeat(t,i) = x(qp.Organize.HeatVented(i)+ qp.Organize.t1States*(t-1));
        end
    end
end
%pull out any dumped cooling
for i = 1:1:length(qp.Organize.CoolVented)
    if qp.Organize.CoolVented(i)>0
        for t = 1:1:n_s
            solution.excessCool(t,i) = x(qp.Organize.CoolVented(i)+ qp.Organize.t1States*(t-1));
        end
    end
end
end%Ends function sort_solution