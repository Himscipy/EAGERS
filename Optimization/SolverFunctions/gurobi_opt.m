function x = gurobi_opt(QP,dt)
global Plant
nS = length(dt);
ic = nnz(QP.Organize.IC);
totalStates = nS*QP.Organize.t1States+ic;
disp = find(QP.Organize.Dispatchable>0);
boolVar = nS*length(disp);
H = diag(QP.H);
f = QP.f;
F = fieldnames(QP.constDemand);
%organize constant demands
constDem = zeros(length(QP.beq),length(disp));
constDem_t = zeros(length(QP.beq),1);
for i = 1:1:length(F)
    if strcmp(F{i},'E')
        net = 'Electrical';
    elseif strcmp(F{i},'H')
        net = 'DistrictHeat';
    elseif strcmp(F{i},'C')
        net = 'DistrictCool';
    end
    for n = 1:1:length(Plant.subNet.(net).nodes) %run through all the nodes in this network
        for t = 1:1:nS
            if QP.constDemand.(F{i}).req((n-1)*nS+t)>0
                for k = 1:1:length(disp)
                    constDem(QP.constDemand.(F{i}).req((n-1)*nS+t),k) = QP.constDemand.(F{i}).load((n-1)*nS+t,disp(k));
                    constDem_t(QP.constDemand.(F{i}).req((n-1)*nS+t),1) = t;
                end
            end
        end
    end
end
%organize which states are dispatchable
dispStates = zeros(totalStates,1);
k = 1;
for i = 1:1:length(disp)
    j = disp(i);
    states = QP.Organize.States{j};
    for t = 1:1:nS
        dispStates(states+(t-1)*QP.Organize.t1States) = k;
        k = k+1;
    end
end

cvx_begin %quiet
% cvx_solver_settings( 'dumpfile', 'test' ) 
% cvx_solver_settings( 'NumericFocus', '3' ) 

    variable x(totalStates)
    variable q(boolVar) binary
    for i = 1:1:totalStates
        C(i,1) = H(i)*x(i)^2 + f(i)*x(i);
    end
    for i = 1:1:totalStates%for all dispatchable gens, add constant cost
        if i<length(disp)
            D(i,1) = QP.constCost(disp(i))*sum(q(nS*(i-1)+1:nS*i));%q from 1 to nS*lengthdisp is boolean related to dispatchable gens,
            %organized by gen1_t to gen1_ns, gen2_t to gen2_ns etc.
        else
            D(i,1) = 0;
        end
    end
    E = sum(C+D);
    
    minimize ( E )
    
    subject to
        QP.A*x <= QP.b;%all inequality constraints
        %add equality constraints
        for r = 1:1:length(QP.beq)
            if constDem_t(r)>0
                QP.Aeq(r,:)*x == QP.beq(r)+ constDem(r,:)*q(constDem_t(r):nS:length(disp)*nS);
            else
                QP.Aeq(r,:)*x == QP.beq(r);
            end
        end
        
        %add upper and lower bounds
        for j = ic+1:1:totalStates
            if dispStates(j)>0
                x(j) <= QP.ub(j)*q(dispStates(j));
                x(j) >= QP.lb(j)*q(dispStates(j));
            else
                x(j) <= QP.ub(j);
                x(j) >= QP.lb(j);
            end
        end
cvx_end