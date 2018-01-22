function sortMC = sortMarginalCost(MarginCost,out,t_max,offGen,dt)
%Converts MarginCost array to 5 columns, sorting by the lowest marginal 
%generation block cost before time t 
%1 = cumulative marginal energy, 
%2 = cumulative marginal cost, used for interpolating to see which blocks are used
%3 = Generator that increases output, 
%4 = time point at which generator increases output, 
%5 = amount that generator increases output (kW)
[nG,nS,n] = size(MarginCost.(out).Capacity.SpinReserve);
MC = zeros(nG*nS*n+1,4);%capacity, cost per kW, generator, time
k = 1;
for i = 1:1:nG
    if isempty(offGen) || i~=offGen
        for t = 1:1:t_max
            if MarginCost.(out).Capacity.SpinReserve(i,t,1)>0
                for j = 1:1:n
                    MC(k+j,1) = MarginCost.(out).Capacity.SpinReserve(i,t,j);
                    MC(k+j,2) = MarginCost.(out).Cost.SpinReserve(i,t,j);
                end
                MC(k+1:k+n,3) = i;
                MC(k+1:k+n,4) = t;
                MC(k+1:k+n,5) = dt(t);
                k = k+n;
            end
        end
    end
end
if k>1
    MC = MC(2:k,:);
    sortMC = zeros(k-1,5);
    [~,I] = sort(MC(:,2));
    sortMC(:,3) = MC(I,3);
    sortMC(:,4) = MC(I,4);
    sortMC(:,5) = MC(I,1);
    sortMC(1,1) = sortMC(1,5)*MC(I(1),5);% (kW*dt) cumulative energy instead of power
    sortMC(1,2) = sortMC(1,5)*MC(I(1),2)*MC(I(1),5);% (kW*$/kW *dt) cumulative cost of the net kJ
    for r = 2:1:k-1
        sortMC(r,1) = sortMC(r-1,1) + sortMC(r,5)*MC(I(r),5);
        sortMC(r,2) = sortMC(r-1,2) + sortMC(r,5)*MC(I(r),2)*MC(I(r),5);%cumulative cost = previous total + size*cost per kW * time
    end
else
    sortMC = [];
end