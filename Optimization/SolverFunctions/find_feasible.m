function [solution,feasible] = find_feasible(gen,qp_0,locked)
%Try removing absorption chiller
n_g = length(gen);
ab_chill = false(1,n_g);
solution.Dispatch = zeros(size(locked));
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Chiller') && strcmp(gen(i).Source,'Heat')
        ab_chill(i) = true;
        locked(:,i) = false;
    end
end
if any(ab_chill)
    locked(:,ab_chill) = false;
    qp = disable_generators(qp_0,locked,[]);%Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution = sort_solution(x,qp);
    end
else
    feasible = 0;
end
if feasible ~=1
    %remove ramping constraints and add them back in until it fails to isolate
    %what constraint is causing it to be infeasible
    %%currently unsure of what the fix is.

    no_ramp = qp_0.Organize.Dispatchable;
    %first check it is feasible without ramping
    QP_1 = remove_ramping(qp_0,no_ramp);
    qp = disable_generators(QP_1,locked,[]);%Disable generators here
    [x,feasible] = call_solver(qp);
    if feasible == 1
        solution = sort_solution(x,qp);
    end
end
% %now try fixing problem
% for i = 1:1:nnz(noRamp)
%     QP_1 = removeRamping(QP_0,noRamp);
%     QP = disable_generators(QP_1,Locked,[]);%Disable generators here
%     [x,Feasible] = call_solver(QP;
%     if Feasible == 1
%         Solution = sort_solution(x,QP);
%     end
% end
end %ends function find_feasible

function qp = remove_ramping(qp,no_ramp)
n_g = length(qp.Organize.Dispatchable);
n_s = length(qp.organize(:,1))-1;
ramp_limit = [];
for i = 1:1:n_g
    if qp.Organize.Dispatchable(i) ==1 && no_ramp(i) ==1 %remove ramping from this generator
        ramp_limit(end+1:end+n_s,1) = (qp.Organize.Ramping(i):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.Ramping(i))';
    end
end
qp.b(ramp_limit) = inf;
qp.b(ramp_limit+1) = inf;
end%Endss function remove_ramping