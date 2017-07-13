function costTerms = GenCosts(Gen,LB,UB,out)
capacity = Gen.Output.Capacity*UB;
efficiency = Gen.Output.(out);
operationRange = find(capacity>=LB);
x = capacity(operationRange);
c = x./efficiency(operationRange); %cost of generator in terms of input
k = find(efficiency(operationRange)==max(efficiency(operationRange)),1,'last');
%% need to add weighting to fitting for uneven spacing of data points

%% FIT A: piecewise convex curve with zero y-intercept (may be same as linear)
costTerms.P = x(k);
costTerms.Convex(1) = c(k)/x(k); %linear segment
if k < length(c)%if index is equal or greater than c, then there is no convex portion
    k2 = max(1,k-1);
    x_i = x(k2:end)-x(k2);%this is (x-P). In the paper alpha is gamma
    y_i = c(k2:end)-costTerms.Convex(1)*x(k2);
    if length(x_i)==2
        x0 = [y_i(2)/x_i(2),0];
    elseif length(x_i)>=3
        x0 = [x_i(2) x_i(2)^2;x_i(3) x_i(3)^2;]\[y_i(2); y_i(3)];
    end
    if length(x_i)>3 || x0(1)<costTerms.Convex(1)
        %convert to constrained linear least squares to quadratic programming problem
        %minimize sum( {d(i) - y(i)}^2 ) where d(i) = c_1*x(i) + c_2*x(i)^2
        %subject to i)  linear slope<= c_1 <= inf  , ii) 0<= c_2 <=inf and iii) c_1*x(end) + c_2*x(end)^2 = y(end)
        %after substitution this becomes: 2*c_1*sum(x(i)*y(i)) + 2*c_2*sum(x(i)^2*y(i)) + c_1^2*sum(x(i)^2) + 2*c_1*c_2*sum(x(i)^3) + c_2^2*sum(x(i)^4) - sum(y(i)^2)
        % put that into quadprog min 0.5*xT*H*x + fT*x form becomes:
        QP.H = 2*[sum(x_i.^2) sum(x_i.^3); sum(x_i.^3) sum(x_i.^4);];
        QP.f = [-2*sum(x_i.*y_i); -2*sum(x_i.^2.*y_i)];
        QP.A = [-1, 0; 0, -1];
        QP.b = [costTerms.Convex(1);0;];
        QP.Aeq = [x_i(end)  x_i(end).^2];
        QP.beq = y_i(end);
        QP.lb = [];
        QP.ub =[];
        QP.solver = 'quadprog';
        [A, cost,Feasible] = callQPsolver(QP);
        if Feasible ==1
            costTerms.Convex(2:3) = A;
        else
            costTerms.Convex = x0;
        end
    else
        costTerms.Convex = x0;
    end
else %if there is no convex portion (Index>= length(c))
    costTerms.Convex(2) = 0;% no quadratic term
    costTerms.Convex(3) = 0; %no quadratic term
end
costTerms.Convex(3) = max(0,costTerms.Convex(3));%ensure a positive H value

%% Fit B: piecewise convex with non-zero y-intercept, first find point I beyond which cost curve is mostly convex
n = length(x);
C = zeros(n,4);
fit = zeros(n,1);
for i = 1:1:n
    I = x(i);% + (i-1)*(x(end) - x(1))/10;
    [C(i,:), fit(i)] = piecewiseQuadratic(I,x,c);
end
%find best fit
[~,index] = min(fit);
costTerms.Intercept = C(index,:);
costTerms.I = x(1) + (index-1)*(x(end) - x(1))/10;
%%Plot
% for i = 1:1:n
%     Y(1:i,i) = x(1:i)*C(i,1)+C(i,4);
%     Y(i+1:n,i) = (x(i+1:n)-x(i))*C(i,2) + (x(i+1:n)-x(i)).^2*C(i,3) + C(i,4)+C(i,1)*x(i);
% end
% figure(1)
% plot(x,c)
% hold on
% plot(x,Y(:,index),'r')



% error = 1;
% best = inf;
% I = (x(end) - x(1))/2 + x(1); %middle of range
% costTerms.I = I;
% while error>1e-2
%     [coef fit] = piecewiseQuadratic(I,x,c);
%     if fit<best
%         costTerms.I = I;
%         costTerms.Intercept = coef;
% 
%         if isinf(best)
%             error = 1;
%         else
%             error = (best-fit)/fit;
%         end
%         best = fit;
%         %take another step in same direction
%     else %backup and take a smaller step
%         
%     end
% end

function [coef, fit] = piecewiseQuadratic(I,x,c)
k = nnz(x<=I);
coef = zeros(1,4);
    
%fit linear segment with c_1*x+c0*x_end
x_i = x(1:k);
y_i = c(1:k);%cost associated with segment
if length(x_i)==1
    coef = [0,0,0,y_i(1)];
else
    x0 = [x_i(1) 1;x_i(2) 1;]\[y_i(1); y_i(2)];
    QP.H = 2*[sum(x_i.^2) sum(x_i); sum(x_i) length(x_i);];
    QP.f = [-2*sum(x_i.*y_i); -2*sum(y_i)];
    QP.Aeq = [x_i(end)  1];
    QP.beq = y_i(end);
    QP.A = [x(end), 1];
    QP.b = c(end);
    QP.lb = [];
    QP.ub = [];
    QP.solver = 'quadprog';
    [A, cost,Feasible] = callQPsolver(QP);
    if Feasible==1
        coef(1) = A(1);
        coef(4) = A(2);
    end
end
fit1 = sum(((x_i*coef(1) + coef(4))- y_i).^2);

x_i = [0; x(k+1:end) - I];
y_i = [0; c(k+1:end) - (I*coef(1) + coef(4))];%cost associated with segment
if k == length(x)
    coef(2) = coef(1);
elseif k == length(x)-1 %not enough points for quadratic segment
    coef(2) = y_i(2)/x_i(2);
else %quadratic segment
    if length(x_i)==2
        x0 = [y_i(2)/x_i(2),0];
    elseif length(x_i)>=3
        x0 = [x_i(2) x_i(2)^2 ;x_i(3) x_i(3)^2;]\[y_i(2); y_i(3)];
    end
    if length(x_i)>3 || x0(1)<coef(1)
        %convert to constrained linear least squares to quadratic programming problem
        %minimize sum( {d(i) - y(i)}^2 ) where d(i) = c_1*x(i) + c2*x(i)^2
        %subject to i)  linear slope<= c_1 <= inf  , ii) 0<= c_2 <=inf and iii) c_1*x(end) + c_2*x(end)^2 = y(end)
        %after substitution this becomes: 2*c_1*sum(x(i)*y(i)) + 2*c_2*sum(x(i)^2*y(i)) + c_1^2*sum(x(i)^2) + 2*c_1*c_2*sum(x(i)^3) + c_2^2*sum(x(i)^4) - sum(y(i)^2)
        % put that into quadprog min 0.5*xT*H*x + fT*x form becomes:
        QP.H = 2*[sum(x_i.^2) sum(x_i.^3); sum(x_i.^3) sum(x_i.^4);];
        QP.f = [-2*sum(x_i.*y_i); -2*sum(x_i.^2.*y_i)];
        QP.Aeq = [x_i(end)  x_i(end).^2];
        QP.beq = y_i(end);
        QP.A = [-1, 0; 0, -1];
        QP.b = [-coef(1);0];
        QP.lb = [];
        QP.ub =[];
        QP.solver = 'quadprog';
        [A,cost, Feasible] = callQPsolver(QP);
        if Feasible==1
            coef(2:3) = A;
        end
    else
        coef(2:3) = x0;
    end
    coef(3) = max(0,coef(3));%ensure positive H value
end
fit = fit1 + sum((coef(2).*x_i + coef(3).*x_i.^2 - y_i).^2);

% UB = max(x);
% Xa = min(x,I);
% Xb = max(x-I,0);
% 
% QP.H = 2*[sum(Xa.^2), sum(Xa.*Xb), sum(Xa.*Xb.^2), sum(Xa); sum(Xa.*Xb), sum(Xb.^2), sum(Xb.^3), sum(Xb); sum(Xa.*Xb.^2), sum(Xb.^3), sum(Xb.^4), sum(Xb.^2); sum(Xa), sum(Xb), sum(Xb.^2), 1;];
% QP.f = -2*[sum(Xa.*c), sum(Xb.*c), sum(Xb.^2.*c), sum(c)];
% QP.Aeq = [I, (UB-I), (UB-I)^2, 1];
% QP.beq = c(end);
% QP.A = [1, -1, 0, 0; 0, -1, 0, 0;];
% QP.b = [0;0;];
% QP.lb = [];
% QP.ub = [];
% QP.solver = 'quadprog';
% [coef, cost,Feasible] = callQPsolver(QP);
% if Feasible==1
%     fit = sum((coef(1)*Xa + coef(2)*Xb + coef(3)*Xb.^2 + coef(4) - c).^2);
% else fit = inf;
% end

% coef(1) = 0;
% X = (x(x>=I)-I)/(x(end)-I);%x normalized to 1
% QP.H = 2*[sum(X.^2), sum(X.^3), sum(X); sum(X.^3), sum(X.^4), sum(X.^2);sum(X), sum(X.^2), 1;];
% QP.f = -2*[sum(X.*c/c(end)), sum(X.^2.*c/c(end)), sum(c/c(end))];%cost normalized to 1
% QP.Aeq = [(1-I/max(x)), (1-I/max(x))^2, 1];
% QP.beq = 1;
% QP.A = [-1, 0, 0; -1, 0, 0;];
% QP.b = [coef(1);0;];
% QP.lb = [];
% QP.ub = [];
% QP.solver = 'quadprog';
% [A, cost,Feasible] = callQPsolver(QP);
% coef = [coef(1), A];
% if Feasible==1
%     fit = sum((coef(1)*min(x,I) + coef(2)*max(x-I) + coef(3)*X.^2 + coef(4) - c).^2);
% else fit = inf;
% end