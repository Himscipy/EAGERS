function K = createCombinations(QP,netDemand,K_C,IC,dt)
global Plant
if isfield(netDemand,'E')
    Out = 'E';
else Out = fieldnames(netDemand); %should only ever be 'C'
    Out = Out{1};
end

nG = length(Plant.Generator);
isCHP = false;
Type = cell(nG,1);
for i = 1:1:nG
    Type(i) = {Plant.Generator(i).Type};
    if strcmp(Type{i},'CHP Generator')
        isCHP = true;
    end
end

if strcmp(Out,'E')
    include = {'Electric Generator';'CHP Generator';'Heater';};%no utilities or storage systems are dispatchable
    if isempty(K_C) && Plant.optimoptions.sequential == 0
        include(end+1) = {'Chiller';};%Need to do something to link absorption chiller to when CHP is on 'Absorption Chiller';
    end
    req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
elseif strcmp(Out,'C')
    include = {'Chiller'};
    req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
elseif strcmp(Out,'H')
    include = {'Heater'};
    req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
end
inc = false(nG,1);
for i = 1:1:nG
    if ismember(Type{i},include)
        states = QP.Organize.States{i};
        if sum(QP.lb(states))>0 %has a lower threshold, and can be shut off
            inc(i) = true;
        end
    end
end
ninc = ~inc;
inc = find(inc);
ninc = find(ninc);
n = length(inc);
K = zeros(2^n,nG); %K is a matrix of all the possible generator on/off combinations 

if ~isempty(ninc)
    for j = 1:1:length(ninc)
        if isempty(K_C) || ~ismember(Type(ninc(j)),'Chiller')
            K(:,ninc(j)) = ninc(j); % all systems that are always included
        end
    end
end

for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=2^n
        K(r+1:r+z,inc(j)) = inc(j);
        r=r+2*z;
    end
end
if ~isempty(K_C)
    nC = length(K_C(:,1));
    K2 = zeros(2^n*nC,nG);
    for i = 1:1:nC
        K2((i-1)*2^n+1:i*2^n,:) = K;
    end
    for j = 1:1:nG
        if ismember(Type(j),'Chiller')
            states = QP.Organize.States{j};
            if sum(QP.lb(states))>0 %has a lower threshold, and can be shut off, so use status determined by K_C
                for i = 1:1:nC
                    K2((i-1)*2^n+1:i*2^n,j) = K_C(i,j);
                end
            end
        end
    end
    K = K2;
end
%% -- %%
%% test if each combination is capable of meeting demand
%if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
%if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
%%could improve this section to check if feasible with network losses (not sure how)
K = UpperLowerLimit(K,QP,req,netDemand,Out,IC,dt);
if isCHP && strcmp(Out,'E') && ~isempty(K)
    req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
    K = UpperLowerLimit(K,QP,req,netDemand,'H',IC,dt);
end
if isempty(K_C) && isfield(Plant.Data,'Demand') && isfield(Plant.Data.Demand,'C') && Plant.optimoptions.sequential == 0 && strcmp(Out,'E') && ~isempty(K)
    req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cooling demand
    K = UpperLowerLimit(K,QP,req,netDemand,'C',IC,dt);
end
%All combinations at this point should be feasible
if isempty(K)
    disp('Zero feasible combinations to test: ERROR')
    K = linspace(1,nG,nG);
end

% Premise (Maybe?), any extra generators beyond what is necessary to meet demand, or
% that are more expensive than grid are unneccessary
%if a combination is feasible, and another combination includes all those
%generators + additional  generators, those other cases can be removed
end%Ends function createCombinations

function K = UpperLowerLimit(K,QP,req,netDemand,Out,IC,dt)
global Plant
nG = length(K(1,:));
limitL = zeros(1,nG);
if strcmp(Out,'H') && Plant.optimoptions.excessHeat
    limitL = limitL - inf;
end
limitU = zeros(1,length(limitL));
for i = 1:1:length(limitL)
    states = QP.Organize.States{i};
    if ~isempty(states) && any(any(QP.Aeq(req,states)~=0))
        I = find(any(QP.Aeq(req,states)~=0));
        lbi = sum(QP.Aeq(req(I),states).*QP.lb(states));
        ubi = sum(QP.Aeq(req(I),states).*QP.ub(states));
        if isfield(Plant.Generator(i).QPform,'Ramp')
            c = ubi/sum(QP.ub(states));
            b = c*Plant.Generator(i).QPform.Ramp.b;
            if isempty(IC) || any(strcmp(Plant.Generator(i).Type,{'Hydro Storage';'Electric Storage';'Thermal Storage';}))%if it is storage
                a = 0;
            else a = c*IC(i);
            end
            if c>0
                lbi = max(lbi,a-b(1)*dt);
                ubi = min(ubi,a+b(2)*dt);
            else %if it consumes this demand (i.e. a chiller when examining electrical demand)
                ubi = max(lbi,a-b(1)*dt);
                lbi = min(ubi,a+b(2)*dt);
            end
        end
        if strcmp(Plant.Generator(i).Type,'Utility')%if it is a utility
            limitU(i) = limitU(i) + ubi;
            limitL(i) = limitL(i) + lbi;
        elseif any(strcmp(Plant.Generator(i).Type,{'Hydro Storage';'Electric Storage';'Thermal Storage';}))%if it is storage, but not a utility
            if ~isempty(IC)
                limitU(i) = limitU(i) + min(IC(i).*QP.Aeq(req,states),ubi);
                chargingSpace = sum((IC(i)-Plant.Generator(i).QPform.Stor.UsableSize).*QP.Aeq(req,states));
                limitL(i) = limitL(i) + max(chargingSpace, lbi);
            end
        else
            limitU(i) = limitU(i) + ubi;
            limitL(i) = limitL(i) + lbi;
        end
    end
end
n = length(K(:,1));
A = (ones(n,1)*limitU);
A(K==0) = 0;
sumUB = sum(A,2);
B = (ones(n,1)*limitL);
B(K==0) = 0;
sumLB = sum(B,2);

% %% Trying to see if we can add transmission losses
% if ~isempty(QP.Organize.Transmission)
%     [sumUB,sumLB] = TransferLoss(QP,A,B,sumUB,sumLB);
% end 

keep = ((sumUB>=netDemand.(Out)) & (sumLB<=netDemand.(Out)));%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
K = K(keep,:);
end%ends function UpperLowerLimit