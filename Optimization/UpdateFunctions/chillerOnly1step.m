function QP = chillerOnly1step(QP,marginal,FirstProfile)
% update the equalities with the correct demand, and scale fuel and electric costs
% EC is the expected end condition at this time stamp (can be empty)
% StorPower is the expected output/input of any energy storage at this timestep (can be empty)
% MinPower and MaxPower define the range of this generator at this timestep
%T is the building temperatures
global Plant
nG = length(Plant.Generator);
abChiller = [];
PartOfDistrictCool = zeros(1,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Chiller')
        PartOfDistrictCool(i) = 1;
        states = QP.Organize.States{i};
        if isfield(Plant.Generator(i).QPform.output,'E')
            QP.f(states) = -marginal.E*Plant.Generator(i).QPform.output.E(:,2);
        elseif isfield(Plant.Generator(i).QPform.output,'H')
            QP.f(states) = -marginal.H*Plant.Generator(i).QPform.output.H(:,2);
            abChiller(end+1) = i;
        end
    elseif strcmp(Plant.Generator(i).Type,'Thermal Storage') && isfield(Plant.Generator(i).QPform.output,'C')
        states = QP.Organize.States{i};
        PartOfDistrictCool(i) = 1;
    end
    %%Do I need to add something to keep the cooling part of a building, or
    %%is that kept automatically?
end

%% set upper limit on ab chiller
if ~isempty(abChiller)
    req = QP.Organize.Balance.DistrictHeat;
    excessHeat = -QP.beq(req);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'CHP Generator';'Heater'}) && FirstProfile(i)>0
            states = QP.Organize.States{i};
            D = 0;
            j = 1;
            if isfield(Plant.Generator(i).QPform,'constDemand')
                excessHeat = excessHeat - Plant.Generator(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
            end
            while D<FirstProfile(i) && j<length(states)
                g = min(FirstProfile(i)-D,QP.ub(states(j)));
                D = D+g;
                excessHeat = excessHeat + g*QP.Aeq(req,states(j));
                j = j+1;
            end
        elseif strcmp(Plant.Generator(i).Type,'Chiller')
            if isfield(Plant.Generator(i).QPform.constDemand,'H')
                QP.constCost(i) = QP.constCost(i) + Plant.Generator(i).QPform.constDemand.H*marginal.H;
            end
            if isfield(Plant.Generator(i).QPform.constDemand,'E')
                QP.constCost(i) = QP.constCost(i) + Plant.Generator(i).QPform.constDemand.E*marginal.E;
            end
        end
    end
    %solve for max abChill output for excess heat
    for k = 1:1:length(abChiller)
        i = abChiller(k);
        states = QP.Organize.States{i};
        for j=1:1:length(states)
            if j == 1
                if excessHeat>(Plant.Generator(i).QPform.constDemand.H+QP.lb(states(1))*(-Plant.Generator(i).QPform.output.H(1,2)))%enough heat to get above LB
                    excessHeat = excessHeat - Plant.Generator(i).QPform.constDemand.H;
                    QP.constCost(i) = 0;
                else
                    excessHeat = 0;
                end
            end
            if QP.ub(states(j))>0
                H = min(excessHeat,QP.ub(states(j))*(-Plant.Generator(i).QPform.output.H(j,2)));
                excessHeat = excessHeat - H;
                freeHeat = H/(-Plant.Generator(i).QPform.output.H(j,2));
                QP.f(states(j)) = (1-freeHeat/QP.ub(states(j)))*QP.f(states(j));%set cost to zero if there is 'free excess heat', but don't reduce the available capacity
            end
        end
    end
end
QP.Organize.HeatVented = [];
%remove all non chillers and cold storage & add cost
QP = rmfield(QP,'constDemand');
QP = disableGenerators(QP,[],PartOfDistrictCool);

%eliminate non district cooling energy balances & other constraints
r = length(QP.b);
xkeep = true(length(QP.f),1);
reqkeep = false(length(QP.beq),1);
rkeep = true(r,1);
reqkeep(QP.Organize.Balance.DistrictCool) = true;
if QP.excessHeat
    xkeep(QP.Organize.HeatVented) = false;
end
%%remove transmission states and constraints
n = length(QP.organize);
nB = length(Plant.Building);
nL = n - nG -nB;
for i = 1:1:nL
    if ~ismember(i,Plant.subNet.DistrictCool.lineNumber)
        xkeep(QP.Organize.States{nG+i}) = false;
        rkeep(QP.Organize.Transmission(i)) = false;
        rkeep(QP.Organize.Transmission(i)+1) = false;
    end
end
if isfield(Plant.subNet,'Electrical') && Plant.optimoptions.SpinReserve
    rkeep(QP.Organize.SpinReserve) = false;
end
for i = 1:1:nB
    rkeep(QP.Organize.Building.r(i)) = false;
end   
QP.H = QP.H(xkeep,xkeep);
QP.f = QP.f(xkeep);
QP.Aeq = QP.Aeq(reqkeep,xkeep);
QP.beq = QP.beq(reqkeep);
if r>0 && nnz(rkeep)>0
    QP.A = QP.A(rkeep,xkeep);
    QP.b = QP.b(rkeep);
else QP.A = [];
    QP.b = [];
end
QP.lb = QP.lb(xkeep);
QP.ub = QP.ub(xkeep);
QP.Organize.t1States = length(QP.f);
QP.Organize.t1Balances = length(Plant.subNet.DistrictCool.nodes);
QP.Organize.Balance = [];
QP.Organize.Balance.DistrictCool = 1;
