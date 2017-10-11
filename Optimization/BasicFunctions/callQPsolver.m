function [x, Cost, Feasible] = callQPsolver(QP,Locked,Enabled,runSetup)
if ~isempty(Locked) || ~isempty(Enabled)
    QP = disableGenerators(QP,Locked,Enabled);%Disable generators here
end
if strcmp(QP.solver,'linprog') && any(any(QP.H - diag(diag(QP.H)))) %not a seperable QP problem
    QP.solver = 'quadprog';
end
if strcmp(QP.solver,'quadprog') && ~license('test','Optimization_Toolbox')
    QP.solver = 'CVX';
end
switch QP.solver
    case 'quadprog'
    %use matlabs linprog or quadprog
    if nnz(QP.H)==0
        options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none'); %,'Display','iter-detailed');% ,'TolFun',1e-10,'TolX',1e-10
        [x,~,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
    else
        options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');%,'TolFun',1e-10,'TolX',1e-10
        [x,~,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options);
    end
    case 'CVX'
        n = length(QP.f);
        if ~isfield(QP,'Organize') || n<24 %calculating generator cost fits
            if ~isempty (QP.lb)
                QP.A = [QP.A;-eye(n);];
                QP.b = [QP.b;QP.lb;];
            end
            if ~isempty (QP.ub)
                QP.ub(isinf(QP.ub)) = 1e3*max(QP.ub(~isinf(QP.ub)));
                QP.A = [QP.A;eye(n);];
                QP.b = [QP.b;QP.ub;];
            end
            cvx_begin quiet
                variable x(n) 
                minimize (0.5*x'*QP.H*x+QP.f'*x)
                subject to
                    QP.Aeq*x == QP.beq;
                    QP.A*x <= QP.b;
            cvx_end
        else
            %% normalize states by ub
            if ~isfield(QP.Organize,'IC')
                ic = 0;
            else ic = max(QP.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
            end
            H = QP.H;
            f = QP.f;
            scale = QP.ub;
            scale(isinf(scale)) = 1e3*max(scale(~isinf(scale)));
            for i=1:1:n
                if QP.ub(i)>0 %for ic = 0 do nothing
                    H(:,i) = QP.H(:,i)*scale(i)^2;
                    f(i) = QP.f(i)*scale(i);
                    if ~isempty(QP.Aeq)
                        QP.Aeq(:,i) = QP.Aeq(:,i)*scale(i);
                    end
                    if ~isempty(QP.A)
                        QP.A(:,i) = QP.A(:,i)*scale(i);
                    end
                    QP.ub(i) = 1;
                end
            end
            %% convert quadprog form to norm
            H = (0.5*diag(H)).^0.5;
            minQcost = min(nonzeros(H));
            b = zeros(n,1);
            for i = 1:1:n
                if QP.f(i)>0
                    if H(i)==0
                        H(i) = 1e-1*minQcost; %a very small quadratic cost so that the linear cost is not NAN
                    end
                    b(i) = f(i)/(-2*H(i));
                end
            end
            H = diag(H);
            %add bound constraints to inequality constraint
            if ~isempty(QP.lb)
                if isempty(QP.A)
                    QP.A = zeros(0,n);
                    QP.b = zeros(0,1);
                end
                I = eye(n);
                I = I(ic+1:end,:); %remove rows that would associate with ic
                QP.A = [QP.A; I;-I];
                QP.b = [QP.b;QP.ub(ic+1:end);-QP.lb(ic+1:end);];
            end
            addpath(genpath('cvx'))
            if runSetup==true %if this is your first time through run the cvx setup
                cvx_setup quiet
            end
            cvx_begin quiet
                variable x(n) nonnegative
                %variable y(n,n) nonnegative
                variable onoff(n) binary
                minimize ((norm((H*y)*x-b)))
                subject to
                    QP.Aeq*x == QP.beq;
                    QP.A*x <= QP.b;
                    QP.lb <= x <= QP.ub;
                    y == zeros(n,n)+diag(onoff);
            cvx_end
            %% convert back to non-normalized
            x = x.*scale;
        end
        if strcmp(cvx_status,'Solved')
            Feasible = 1;
        else Feasible =0;
        end
    case 'SPDT3'
    
    case 'linprog'
        %piecewise fit of any quadratic terms
        % error relative to quadprog depends upon the # of linear segments
        n = 5; %doubling the segments = 1/2 the error
        Hd = diag(QP.H);
        index = nonzeros((1:length(Hd))'.*(Hd>0));
        xL = length(Hd)+(n-1)*length(index); %new length of x vector;
        f = zeros(xL,1);
        Aeq = zeros(length(QP.beq),xL);
        beq = QP.beq;
        A = zeros(length(QP.b),xL);
        b = QP.b;
        lb = zeros(xL,1);
        ub = zeros(xL,1);
        k = 0;
        for i = 1:1:length(Hd)
            if any(index==i) %split state into n linear segments
                r = QP.ub(i) - QP.lb(i);
                x = linspace(QP.lb(i)+r/n,QP.ub(i),n)';
                f_i = QP.f(i) + QP.H(i,i)*x;%linear slope of each quadratic section
                for j = 1:1:n
                    k = k+1;
                    f(k) = f_i(j);
                    Aeq(:,k) = QP.Aeq(:,i);
                    A(:,k) = QP.A(:,i);
                    if j ==1
                        lb(k) = QP.lb(i);
                        ub(k) = QP.lb(i)+r/n;
                    else
                        lb(k) = 0;
                        ub(k) = r/n;
                    end
                end
            else %copy state from before
                k = k+1;
                f(k) = QP.f(i);
                Aeq(:,k) = QP.Aeq(:,i);
                A(:,k) = QP.A(:,i);
                lb(k) = QP.lb(i);
                ub(k) = QP.ub(i);
            end
        end
        options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none');
        [x2,~,Feasible] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        if Feasible == 1 %recover x back into original length
            x = zeros(length(Hd),1);
            k = 0;
            for i = 1:1:length(x)
                if any(index==i) %sum the split state
                    x(i) = sum(x2(k+1:k+n));
                    k = k+n;
                else
                    x(i) = x2(k+1);
                    k = k+1;
                end
            end
        end
    case 'PredictorCorrector'
        n = length(QP.f);
        %% normalize states by ub
        if ~isfield(QP.Organize,'IC')
            ic = 0;
        else ic = max(QP.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
        end
%         scale = QP.ub;
%         scale(isinf(scale)) = 1e3*max(scale(~isinf(scale)));
%         for i=1:1:n
%             if QP.ub(i)>0 %for ic = 0 do nothing
%                 QP.H(:,i) = QP.H(:,i)*scale(i)^2;
%                 QP.f(i) = QP.f(i)*scale(i);
%                 if ~isempty(QP.Aeq)
%                     QP.Aeq(:,i) = QP.Aeq(:,i)*scale(i);
%                 end
%                 if ~isempty(QP.A)
%                     QP.A(:,i) = QP.A(:,i)*scale(i);
%                 end
%                 QP.ub(i) = 1;
%             end
%         end
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
%         x = x.*scale;
end
if Feasible ==1
    Cost = 0.5*x'*QP.H*x + x'*QP.f;
else Cost = 0;
end
if isfield(QP,'organize')
    [m,n] = size(QP.organize);
    nG = length(QP.constCost);
    nS = m-1;
    GenDisp = zeros(m,n);
    if Feasible ==1
        %%Getting data for analysis
        for i = 1:1:nG
            if QP.Organize.Hydro(i) == 1
                for t = 1:1:nS+1
                    %Get SOC of each generator into a matrix for all time steps
                    SOCGen(t,i) = x(QP.organize{t,i}+1,1);
                    global Plant
                    Plant.SOCSolvedQP = SOCGen;
                end
            end
        end
        for i = 1:1:n
            if i > nG && i <= n-nG
                for t = 1:1:nS+1
                    if ~isempty(QP.organize{t,i})
                        Trans(t,i-nG) = sum(x(QP.organize{t,i}+1)); %down (positive) lines
                        Trans(t,i) = sum(x(QP.organize{t,i}+2)); %up (negative) lines
                    end
                    if t == nS+1
                        Plant.Trans = Trans;
                    end
                end
            end
        end
        %end analysis data
        for i = 1:1:n
            if i<=nG && isfield(QP,'Renewable') && any(QP.Renewable(:,i)~=0)
                GenDisp(1,i) = RenewableOutput(i,[],'Actual');
                GenDisp(2:end,i) = QP.Renewable(:,i);
            else
                for t = 1:1:nS+1
                    if ~isempty(QP.organize{t,i})
                        GenDisp(t,i) = sum(x(QP.organize{t,i})); %record this combination of outputs (SOC for storage)
                    end
                end
            end
        end
        GenDisp(abs(GenDisp)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
    end
    x = GenDisp;
end
end%Ends callQPsolver
