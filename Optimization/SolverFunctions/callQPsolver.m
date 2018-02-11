function [x,Feasible] = callQPsolver(QP)
if strcmp(QP.solver,'linprog') && any(any(QP.H - diag(diag(QP.H)))) %not a seperable QP problem
    QP.solver = 'quadprog';
end
if strcmp(QP.solver,'quadprog') && ~license('test','Optimization_Toolbox')
%     QP.solver = 'Gurobi';
    QP.solver = 'qpip';
end
if isempty(QP.f)
    x = [];
    Feasible = 1;
else
    switch QP.solver
        case 'quadprog'
        %use matlabs linprog or quadprog
        if ~any(QP.H)
            options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none'); %,'Display','iter-detailed');% ,'TolFun',1e-10,'TolX',1e-10
            [x,~,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
        else
            options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');%,'TolFun',1e-10,'TolX',1e-10
            [x,~,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options);
            if Feasible ==0
                Feasible = 1;
%                 disp('Max iterations exceeded')
            elseif Feasible == -2 || Feasible == 3
                if ~isempty(x) && isfield(QP.Organize.Balance,'DistrictHeat')%change scale of heat vented upper bound
                    [m,n] = size(QP.organize);
                    nS = max(1,m-1);
                    req = nonzeros((1:1:length(QP.Aeq(:,1)))'.*(QP.Aeq(:,QP.Organize.HeatVented)~=0));
                    req = (req:QP.Organize.t1Balances:(req+(nS-1)*QP.Organize.t1Balances))';
                    error = QP.Aeq*x-QP.beq;
                    heatventStates = (QP.Organize.HeatVented:QP.Organize.t1States:(QP.Organize.HeatVented+(nS-1)*QP.Organize.t1States))';
                    QP.ub(heatventStates) = max(10,x(heatventStates) + 2*error(req));
                    [x,~,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options);
                    %error = QP.Aeq*x-QP.beq; %% can use this error to determine what energy balance is not met. can help in filter gen
%                     if Feasible ==1
%                         disp('Converged, but heat state needed rescaling')
%                     end
                end
            end
        end
        case 'Gurobi'
            [x,Feasible] = gurobi_opt(QP);
        case 'qpip'
            %move upper and lower bounds into inequality constraints
            n = length(QP.f);
            if ~isfield(QP,'Organize')|| ~isfield(QP.Organize,'IC')
                ic = 0;
            else ic = max(QP.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
            end
            if ~isempty(QP.lb)
                if isempty(QP.A)
                    QP.A = zeros(0,n);
                    QP.b = zeros(0,1);
                end
                QP.ub(isinf(QP.ub)) = 1e3*max(QP.ub(~isinf(QP.ub)));%% avoid inf bounds!!
                I = [zeros(n-ic,ic),eye(n-ic)];
                QP.A = [QP.A; I;-I];
                QP.b = [QP.b;QP.ub(ic+1:end);-QP.lb(ic+1:end);];
            end
            [x,Feasible,stats] = qpip(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq);
            x = x.x;
        case 'SPDT3'
            %unclear how to convert into the sparse conic problem definition
        case 'MIPS'
            %see http://www.pserc.cornell.edu/matpower/docs/MIPS-manual-1.2.2.pdf
            %and https://github.com/MATPOWER/mips
        case 'PredictorCorrector'
            n = length(QP.f);
            %% normalize states by ub
            if ~isfield(QP.Organize,'IC')
                ic = 0;
            else ic = max(QP.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
            end
            scale = QP.ub;
            scale(isinf(scale)) = 1e3*max(scale(~isinf(scale)));
            for i=1:1:n
                if QP.ub(i)>0 %for ic = 0 do nothing
                    QP.H(:,i) = QP.H(:,i)*scale(i)^2;
                    QP.f(i) = QP.f(i)*scale(i);
                    if ~isempty(QP.Aeq)
                        QP.Aeq(:,i) = QP.Aeq(:,i)*scale(i);
                    end
                    if ~isempty(QP.A)
                        QP.A(:,i) = QP.A(:,i)*scale(i);
                    end
                    QP.ub(i) = 1;
                end
            end
            %add bound constraints to inequality constraint
            if ~isempty(QP.lb)
                if isempty(QP.A)
                    QP.A = zeros(0,n);
                    QP.b = zeros(0,1);
                end
                QP.ub(isinf(QP.ub)) = 1e3*max(QP.ub(~isinf(QP.ub)));%% avoid inf bounds!!
                I = [zeros(n-ic,ic),eye(n-ic)];
                QP.A = [QP.A; I;-I];
                QP.b = [QP.b;QP.ub(ic+1:end);-QP.lb(ic+1:end);];
            end
            [r,~] = size(QP.A);
            [req,~] = size(QP.Aeq);
            % initial values for states, langragian coef, and slack variables
            x = zeros(n,1);
            y = ones(req,1);
            z = ones(r,1);
            s = ones(r,1);
            [x,iterations,Feasible] = pcQPgen(QP.H,QP.f,QP.A',QP.b,QP.Aeq',QP.beq,x,y,z,s);
            %% convert back to non-normalized
            x = x.*scale;
    end
end
end%Ends callQPsolver
