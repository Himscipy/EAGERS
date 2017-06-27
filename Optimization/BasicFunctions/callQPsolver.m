function [x, Cost, Feasible] = callQPsolver(QP)
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
        options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none');
        [x,~,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
    else
        options2 = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
        [x,~,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2);
    end
    case 'CVX'
        n = length(QP.f);
        if ~isfield(QP,'Organize') || n<24 %calculating generaor cost fits
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
            cvx_begin quiet
                variable x(n) nonnegative
                minimize ((norm(H*x-b)))
                subject to
                    QP.Aeq*x == QP.beq;
                    QP.A*x <= QP.b;
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
end
if Feasible ==1
    Cost = 0.5*x'*QP.H*x + x'*QP.f;
else Cost = 0;
end
