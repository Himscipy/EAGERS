function [x,feasible] = call_solver(qp)
if strcmp(qp.solver,'ANN')
    qp.solver = 'quadprog';
end
if strcmp(qp.solver,'linprog') && any(any(qp.H - diag(diag(qp.H)))) %not a seperable QP problem
    qp.solver = 'quadprog';
end
if strcmp(qp.solver,'quadprog') && ~license('test','Optimization_Toolbox')
    qp.solver = 'qpip';
end
if isempty(qp.f)
    x = [];
    feasible = 1;
else
    switch qp.solver
        case 'quadprog'
        %use matlabs linprog or quadprog
        if ~any(qp.H)
            options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none'); %,'Display','iter-detailed');% ,'TolFun',1e-10,'TolX',1e-10
            [x,~,feasible] = linprog(qp.f,qp.A,qp.b,qp.Aeq,qp.beq,qp.lb,qp.ub,[],options); 
        else
            options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');%,'TolFun',1e-10,'TolX',1e-10
            [x,~,feasible] = quadprog(qp.H,qp.f,qp.A,qp.b,qp.Aeq,qp.beq,qp.lb,qp.ub,[],options);
            if feasible ==0
                feasible = 1;
%                 disp('Max iterations exceeded')
            elseif feasible == -2 || feasible == 3
                if ~isempty(x) && isfield(qp.Organize.Balance,'DistrictHeat')%change scale of heat vented upper bound
                    [m,n] = size(qp.organize);
                    nS = max(1,m-1);
                    req = nonzeros((1:1:length(qp.Aeq(:,1)))'.*(qp.Aeq(:,qp.Organize.HeatVented)~=0));
                    req = (req:qp.Organize.t1Balances:(req+(nS-1)*qp.Organize.t1Balances))';
                    error = qp.Aeq*x-qp.beq;
                    heatventStates = (qp.Organize.HeatVented:qp.Organize.t1States:(qp.Organize.HeatVented+(nS-1)*qp.Organize.t1States))';
                    qp.ub(heatventStates) = max(10,x(heatventStates) + 2*error(req));
                    [x,~,feasible] = quadprog(qp.H,qp.f,qp.A,qp.b,qp.Aeq,qp.beq,qp.lb,qp.ub,[],options);
                    %error = QP.Aeq*x-QP.beq; %% can use this error to determine what energy balance is not met. can help in filter gen
%                     if Feasible ==1
%                         disp('Converged, but heat state needed rescaling')
%                     end
                end
            end
        end
        case 'Gurobi'
            [x,feasible] = gurobi_opt(qp);
        case 'qpip'
            [x,feasible] = qpip(qp);
        case 'SPDT3'
            %unclear how to convert into the sparse conic problem definition
        case 'MIPS'
            %see http://www.pserc.cornell.edu/matpower/docs/MIPS-manual-1.2.2.pdf
            %and https://github.com/MATPOWER/mips
        case 'PredictorCorrector'
            n = length(qp.f);
            %% normalize states by ub
            if ~isfield(qp.Organize,'IC')
                ic = 0;
            else
                ic = max(qp.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
            end
            scale = qp.ub;
            scale(isinf(scale)) = 1e3*max(scale(~isinf(scale)));
            for i=1:1:n
                if qp.ub(i)>0 %for ic = 0 do nothing
                    qp.H(:,i) = qp.H(:,i)*scale(i)^2;
                    qp.f(i) = qp.f(i)*scale(i);
                    if ~isempty(qp.Aeq)
                        qp.Aeq(:,i) = qp.Aeq(:,i)*scale(i);
                    end
                    if ~isempty(qp.A)
                        qp.A(:,i) = qp.A(:,i)*scale(i);
                    end
                    qp.ub(i) = 1;
                end
            end
            %add bound constraints to inequality constraint
            if ~isempty(qp.lb)
                if isempty(qp.A)
                    qp.A = zeros(0,n);
                    qp.b = zeros(0,1);
                end
                qp.ub(isinf(qp.ub)) = 1e3*max(qp.ub(~isinf(qp.ub)));%% avoid inf bounds!!
                I = [zeros(n-ic,ic),eye(n-ic)];
                qp.A = [qp.A; I;-I];
                qp.b = [qp.b;qp.ub(ic+1:end);-qp.lb(ic+1:end);];
            end
            [r,~] = size(qp.A);
            [req,~] = size(qp.Aeq);
            % initial values for states, langragian coef, and slack variables
            x = zeros(n,1);
            y = ones(req,1);
            z = ones(r,1);
            s = ones(r,1);
            [x,iterations,feasible] = pcQPgen(qp.H,qp.f,qp.A',qp.b,qp.Aeq',qp.beq,x,y,z,s);
            %% convert back to non-normalized
            x = x.*scale;
    end
end
end%Ends call_solver