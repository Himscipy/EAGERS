function [x,k,feas] = pcQPgen(G,g,C,d,A,b,x,y,z,s)
% This function solves a QP problem of the form
% min 0.5x’Gx+g’x
% s.t. A’x=b
% C’x>=d
% using the Predictor?Corrector (PC) method.
%
% Input :
% x:starting point for the vector x (nx1 vector)
% y:starting point for the vector y (nAx1 vector) where nA = # of equality constraints
% z:starting point for the vector z (nCx1 vector) where nC = # of inequality constraints
% s:starting point for the slack vector s (nCx1 vector)
% G:Hessian (nxn matrix)
% g:(nx1 vector)
% C: left hand side of inequality constraints(nxnC matrix)
% d:right hand side of inequality constraints (nCx1 vector)
% A:left hand side of equality constraints (nxnA matrix)
% b:right hand side of equality constraints (nAx1 vector)
% where nC and nA are the numbers of inequality and equality constraints.

% Output:
% xstop :solution x
% k : number of iterations used
%
% Thomas Res low Kr¨uth , s 0 2 1 8 9 8
%dampening factor eta
eta = 0.95;
%residuals are computed
[mA,nA] = size(A);
[mC,nC] = size(C);
e = ones(nC,1);
rL = G*x+g-A*y-C*z;
rA = -A'*x+b;
rC = -C'*x+s+d;
rsz = s.*z;
mu = sum(z.*s)/nC;
%k: number of iterations, epsilon: tolerances
k = 0;
maxk = 200;
epsL = 1e-10;
epsA=1e-10;
epsC=1e-10;
epsmu=1e-10;
con = true; %continue
while con
    %Solve system with a Newton like method/Factorizarion
    lhs = [G,-A,-C;-A',spalloc(nA,nA,0),spalloc(nA,nC,0);-C',spalloc(nC,nA,0),sparse(-diag(s./z))];
    [L,D,P]=ldl(lhs);
    rhs = [-rL;-rA;-rC+rsz./z];
    dxyz_a = P*(L'\(D\(L\(P'*rhs))));
    dx_a = dxyz_a(1:length(x));
    dy_a = dxyz_a(length(x)+1:length(x)+length(y));
    dz_a = dxyz_a(length(x)+length(y)+1:length(x)+length(y)+length(z));
    ds_a = -((rsz+s.*dz_a)./z);
    %Compute alphaaff
    alpha_a = 1;
    idx_z = find(dz_a<0);
    if (isempty(idx_z)==0)
        alpha_a = min(alpha_a,min(-z(idx_z)./dz_a(idx_z)));
    end
    idx_s = find(ds_a<0);
    if(isempty(idx_s)==0)
        alpha_a= min(alpha_a,min(-s(idx_s)./ds_a(idx_s)));
    end
    %Compute the affine duality gap
    mua = ((z+alpha_a*dz_a)'*(s+alpha_a*ds_a))/nC;
    %Compute the centering parameter
    sigma = (mua/mu)^3;
    %Solve system
    rsz = rsz + ds_a.*dz_a - sigma*mu*e;
    rhs = [-rL;-rA;-rC+rsz./z];
    dxyz = P*(L'\(D\(L\(P'*rhs))));
    dx = dxyz(1:length(x));
    dy = dxyz(length(x)+1:length(x)+length(y));
    dz = dxyz(length(x)+length(y)+1:length(x)+length(y)+length(z));
    ds = -((rsz+s.*dz)./z);
    %Compute alpha
    alpha = 1;
    idx_z = find(dz<0);
    if (isempty(idx_z)==0)
        alpha = min(alpha,min(-z(idx_z)./dz(idx_z)));
    end
    idx_s = find(ds<0);
    if (isempty(idx_s)==0)
        alpha = min(alpha,min(-s(idx_s)./ds(idx_s)));
    end
    %Update x, z, s
    x = x + eta*alpha*dx;
    y = y + eta*alpha*dy;
    z = z + eta*alpha*dz;
    s= s + eta*alpha*ds;
    
    %Update rhs
    rL = G*x+g-A*y-C*z;
    rA = -A'*x+b;
    rC = -C'*x+s+d;
    rsz = s.*z;
    mu = sum(z.*s)/nC;
    
    
    if k==maxk
        con = false;
        feas = 0;
    else
        k = k+1;
        if norm(rL)<epsL && norm(rA)<epsA && norm(rC)<epsC && abs(mu)<epsmu
            con = false;
            feas = 1;
        end
    end  
end
