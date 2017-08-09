function K = createCombinations(QP,netDemand)
global Plant
if~isempty(Plant.Building)
    %% This is temporary to only do electric dispatch
    Outs = {'E'};
else
    Outs = fieldnames(Plant.Data.Demand);
end
nG = length(Plant.Generator);
Type = cell(nG,1);
for i = 1:1:nG
    Type(i) = {Plant.Generator(i).Type};
end

isCHP = false;
Outs = Outs(~strcmp('W',Outs)); %no dams are dispatchable
if any(strcmp('E',Outs)) && any(strcmp('H',Outs))%check if there is a CHP system
    for i = 1:1:nG
        if strcmp(Type{i},'CHP Generator')
            isCHP = true;
            break
        end
    end
    if isCHP
        Outs = Outs(~strcmp('H',Outs)); %heaters are part of electricity case
    end
end
if any(strcmp('E',Outs)) && any(strcmp('C',Outs)) && Plant.optimoptions.sequential == 0
    Outs = Outs(~strcmp('C',Outs)); %Chillers are part of electricity case
end
K = zeros(0,nG);
lines = 0;
for s = 1:1:length(Outs)
    if strcmp(Outs{s},'E')
        include = {'Electric Generator';'CHP Generator';};%no utilities are dispatchable
        if Plant.optimoptions.sequential == 0
            include(end+1) = {'Chiller';};%no utilities are dispatchable ___ Need to do something to link absorption chiller to when CHP is on 'Absorption Chiller';
        end
        if isCHP
            include(end+1) = {'Heater'};
        end
        req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    elseif strcmp(Outs{s},'C')
        include = {'Chiller'};
        req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    elseif strcmp(Outs{s},'H')
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
    K = [K; zeros(2^n,nG)]; %K is a matrix of all the possible generator on/off combinations 
    if ~isempty(inc) && ~isempty(ninc)
        K(lines+1:end,ninc) = ones(2^n,1)*ninc'; % all systems that are always included
    end
    for j = 1:1:n %all combinations of generators are listed below
        z = 2^(n-j);
        r=0;
        while r+z<=2^n
            K(lines+r+1:lines+r+z,inc(j)) = inc(j);
            r=r+2*z;
        end
    end
    
    %% -- %%
    %% test if each combination is capable of meeting demand
    %if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
    %if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
    %%could improve this section to check if feasible with network losses (not sure how)
    limitL = zeros(1,nG);
    if strcmp(Outs{s},'H') && Plant.optimoptions.excessHeat
        limitL = limitL - inf;
    end
    K = UpperLowerLimit(K,lines,QP,req,netDemand.(Outs{s}),limitL);
    if isCHP && strcmp(Outs{s},'E') && ~isempty(K)
        req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
        limitL = zeros(1,nG);
        if Plant.optimoptions.excessHeat
            limitL = limitL - inf;
        end
        K = UpperLowerLimit(K,lines,QP,req,netDemand.H,limitL);
    end
    if isfield(Plant.Data,'Demand') && isfield(Plant.Data.Demand,'C') && Plant.optimoptions.sequential == 0 && strcmp(Outs{s},'E') && ~isempty(K)
        req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cooling demand
        limitL = zeros(1,nG);
        K = UpperLowerLimit(K,lines,QP,req,netDemand.C,limitL);
    end
    lines = length(K(:,1));
end
end%Ends function createCombinations

function K = UpperLowerLimit(K,lines,QP,req,netDemand,limitL)
limitU = zeros(1,length(limitL));
for i = 1:1:length(limitL)
    states = QP.Organize.States{i};
    if~isempty(states)
        if any(isinf(QP.ub(states)))
            if any(QP.Aeq(req,states)~=0)
                limitU(i) = inf;
            end
        else
            limitU(i) = limitU(i) + sum(max(0,QP.Aeq(req,states)).*QP.ub(states));
            limitU(i) = limitU(i) - sum(min(0,QP.Aeq(req,states)).*QP.ub(states));%add anything thats an electric load (i.e. chiller))
        end
        if any(isinf(QP.lb(states)))
            if any(QP.Aeq(req,states)~=0)
                limitL(i) = -inf;
            end
        else
            limitL(i) = limitL(i) + sum(max(0,QP.Aeq(req,states).*QP.lb(states)));
            limitL(i) = limitL(i) - sum(min(0,QP.Aeq(req,states).*max(0,QP.lb(states))));
        end
    end
end
n = length(K(:,1))-lines;
A = (ones(n,1)*limitU);
A(K(lines+1:end,:)==0) = 0;
sumUB = sum(A,2);
B = (ones(n,1)*limitL);
B(K(lines+1:end,:)==0) = 0;
sumLB = sum(B,2);
keep = [true(lines,1); ((sumUB>=netDemand) & (sumLB<=netDemand))];%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
K = K(keep,:);
end%ends function UpperLowerLimit
