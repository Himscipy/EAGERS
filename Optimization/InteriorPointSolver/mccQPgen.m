function[xstop,ystop,zstop,sstop,k] = mccQPgen(x,y,z,s,G,g,C,d,A,b)
%??????????????????????????????????????????????????????????
% mccQPgen.m
%
%This function solves a QP problem of the form
% min 0.5x’Gx+g’x
%s.t. A’x = b
% C’x>=d
% using the Multiple Centrality Corrections (MCC) method .
%
% Input :
% x:starting point for the vector x (nx1 vector)
% y:starting point for the vector x (nAx1 vector)
% z:starting point for the vector x (nCx1 vector)
% s:starting point for the slack ?vector x (nCx1 vector)
% G:Hessian (nxn matrix)
% g:(nx1 vector)
% C: left hand side of in equality constraints(nxnC matrix)
% d:right hand side of inequality constraints (nCx1 vector)
% A:left hand side of equality constraints (nxnA matrix)
% b:right hand side of equality constraints (nAx1 vector)
% where nC and nA are the numbers of inequality and equality
% constraints .
% Output:
% xstop :solution x
% ystop
% zstop
% sstop
% k : number of iterations used
%
% Thomas Res low Kr¨uth , s 0 2 1 8 9 8
%??????????????????????????????????????????????????????????
Bmin = 0.1; Bmax = 10; dalpha = 0.1; gamma = 0.1;
eta = 0.99995;
%residuals are computed
[mA,nA] = size(A);
[mC,nC] = size(C);
e = ones(nC,1);
rL = G*x + g - A*y - C*z ;
rA = -A'*x + b;
rC = -C'*x + s+d;
rsz = s.*z;
mu = sum(z.*s)/nC;
%k: number of iterations, epsilon: tolerances
k = 0 ;
maxk = 200;
epsL = 1e-10;epsA = 1e-10; epsC = 1e-10; epsmu = 1e-10;
while(k<=maxk && norm(rL)>=epsL && norm(rA)>=epsA && norm(rC)>=epsC && abs(mu)>=epsmu)
    %Solve system with a Newton?like method/Factorizarion
    lhs = [G,-A,-C;-A',spalloc(nA,nA,0),spalloc(nA,nC,0);-C',spalloc(nC,nA,0),sparse(-diag(s./z))];
    [L,D,P] = ldl(lhs);
    rhs = [-rL;-rA;-rC+rsz./z];
    dxyza = P*(L’\(D\(L\(P'*rhs))));
    dza = dxyza(length(x)+length(y)+1:length(x)+length(y)+length(z));
    dsa = -((rsz+s.*dza)./z);
    %Compute alpha aff
    alphaa = 1;
    idxz = find(dza<0);
    if (isempty(idxz)==0)
        alphaa = min(alphaa,min(-z(idxz)./dza(idxz)));
    end
    idxs = find(dsa<0);
    if (isempty(idxs)==0)
        alphaa = min(alphaa,min(-s(idxs)./dsa(idxs)));
    end
    %Compute the affine duality gap
    mua = ((z+alphaa*dza)'*(s+alphaa*dsa))/nC;
    %Compute the centering parameter
    sigma = (mua/mu)^3;
    mut = sigma*mu;
    rsz = rsz+dsa.*dza-mut*e;
    rhs = [-rL;-rA;-rC+rsz./z];
    dxyzp = P*(L'\(D\(L\(P'*rhs))));
    dxp = dxyzp(1:length(x));
    dyp = dxyzp(length(x)+1:length(x)+length(y));
    dzp = dxyzp(length(x)+length(y)+1:length(x)+length(y)+length(z));
    dsp = -((rsz+s.*dzp)./z);
    %Computealpha
    alphap = 1;
    idxz = find(dzp<0);
    if (isempty(idxz)==0)
        alphap = min(alphap,min(-z(idxz)./dzp(idxz)));
    end
    idxs = find(dsp<0);
    if (isempty(idxs)==0)
        alphap=min(alphap,min(-s(idxs)./dsp(idxs)));
    end
    %Modifiedcentering directions
    alphah = min(alphap+dalpha,1);
    %Maximum number of corrections
    K=2;
    %Corrections
    j=1;
    while(j<=K)
        %Compute trial point
        zt = z+alphah*dzp;
        st = s+alphah*dsp;
        %Define target
        vthilde = st.*zt;
        vtmp= vthilde;
        for i =1:length(vtmp)
            if (vtmp(i)<Bmin*mut)
                vtmp(i) = Bmin*mut;
            elseif (vtmp(i)>Bmax*mut)
                vtmp(i) = Bmax*mut;
            end
        end
        vt = vtmp;
        %Compute corrector
        rsz = -(vt-vthilde);
        for i=1:length(rsz)
            if (rsz(i)<(-Bmax*mut))
                rsz(i) = -Bmax*mut;
            end
        end
        rhs =[zeros(length(rL),1);zeros(length(rA),1);rsz./z];
        dxyzm = P*(L'\(D\(L\(P'*rhs))));
        dxm = dxyzm(1:length(x));
        dym = dxyzm(length(x)+1:length(x)+length(y));
        dzm = dxyzm(length(x)+length(y)+1:length(x)+length(y)+length(z));
        dsm = -((rsz+s.*dzm)./z);
        %Composite direction
        dx = dxp+dxm;
        ds = dsp+dsm;
        dy = dyp+dym;
        dz = dzp+dzm;
        %Compute alpha
        alpha = 1;
        idxz = find(dz<0);
        if (isempty(idxz)==0)
            alpha = min(alpha,min(-z(idxz)./dz(idxz)));
        end
        idxs = find(ds<0);
        if (isempty(idxs)==0)
            alpha = min(alpha,min(-s(idxs)./ds(idxs)));
        end
        %Test for improvement
        if (alpha>=alphap+gamma*dalpha)
            j=j+1;
            dxp = dx;
            dsp = ds;
            dzp = dz;
            dyp = dy;
            alphap = alpha;
            alphah = min(alphap+dalpha,1);
        else
            dx = dxp;
            ds = dsp;
            dz = dzp;
            dy = dyp;
            j=10;
        end
    end %inner while?loop
    %Update x,y,z,s
    x = x+eta*alpha*dx;
    y = y+eta*alpha*dy;
    z = z+eta*alpha*dz;
    s = s+eta*alpha*ds;
    k = k+1;
    %Update rhs
    rL = G*x+g-A*y-C*z;
    rA = -A'*x+b;
    rC = -C'*x+s+d;
    rsz = s.*z;
    mu = sum(z.*s)/nC;
end
%Output
xstop = x;
ystop =y;
zstop = z;
sstop = s;
