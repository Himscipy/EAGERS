function QP = chillerOnly1step(QP,marginal)
% update the equalities with the correct demand, and scale fuel and electric costs
% EC is the expected end condition at this time stamp (can be empty)
% StorPower is the expected output/input of any energy storage at this timestep (can be empty)
% MinPower and MaxPower define the range of this generator at this timestep
%T is the building temperatures
global Plant
nG = length(Plant.Generator);

PartOfDistrictCool = zeros(1,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Chiller') || (strcmp(Plant.Generator(i).Type,'Thermal Storage') && isfield(Plant.Generator(i).QPform.output,'C'))
        PartOfDistrictCool(i) = 1;
        states = QP.Organize.States{i};
        if isfield(Plant.Generator(i).QPform.output,'E')
            QP.f(states) = -marginal.Electrical*Plant.Generator(i).QPform.output.E(1,1);
        elseif isfield(Plant.Generator(i).QPform.output,'H')
            QP.f(states) = -marginal.DistrictHeat*Plant.Generator(i).QPform.output.H(1,1);
        end
    end
    %%Do I need to add something to keep the cooling part of a building, or
    %%is that kept automatically
end
%remove all non chillers and cold storage & add cost
QP = rmfield(QP,'constDemand');
QP = disableGenerators(QP,[],PartOfDistrictCool);

for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'constDemand')
        if isfield(Plant.Generator(i).QPform.constDemand,'H')
            QP.constCost(i) = QP.constCost(i) + Plant.Generator(i).QPform.constDemand.H*marginal.DistrictHeat;
        end
        if isfield(Plant.Generator(i).QPform.constDemand,'E')
            QP.constCost(i) = QP.constCost(i) + Plant.Generator(i).QPform.constDemand.E*marginal.Electrical;
        end
    end
end

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
