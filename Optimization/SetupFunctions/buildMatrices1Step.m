function QP = buildMatrices1Step
%builds constant matrices for step-by-step optimization
%similar to multi-step except:
%energy storage looks like a generator with 1 state (can be pos or neg)
%no ramping constraints
%Fit A includes energy storage and uses the fit with zero y-intercept
%Fit B does not include energy storage and uses non-zero intercept
%Demands, upper/lower bounds, and utility costs must updated prior to optimization
global Plant 
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.(networkNames{net}).lineNames);
end
nL = sum(nLinet);

Organize.States = cell(1,nG+nL+nB);
Organize.Dispatchable = zeros(1,nG);
Organize.Transmission = zeros(nL,1);
Organize.Hydro = false(1,nG);
Organize.Building.r = zeros(nB,1);
Organize.Building.req = zeros(nB,1);
Organize.Out_vs_State = cell(1,nG+nL+nB);
QP.organize = cell(1,nG+nL+nB);
QP.constCost = zeros(1,nG);
[QP,Organize,lb,ub] = countStates1Step(QP,Organize,nG,nB,nL);
H = zeros(Organize.t1States,1);
QP.f = zeros(Organize.t1States,1);
QP.ub = ub;
QP.lb = lb;

Organize = equaltiyConstraints1Step(Organize);
QP.Aeq = zeros(Organize.t1Balances,Organize.t1States);
QP.beq = zeros(Organize.t1Balances,1);

Organize = inequalityConstraints1Step(Organize,nG,nB);
QP.A = zeros(Organize.t1ineq,Organize.t1States);
QP.b = zeros(Organize.t1ineq,1);

%%Call functions to fill in these matrices
[QP.Aeq,QP.beq,H,QP.f] = EqualityConstraints(QP.Aeq,QP.beq,H,QP.f,Organize);
QP.Aeq = LineEqualities(QP.Aeq,Organize);
QP.A = StorageInequalities(QP.A,Organize);
QP.A = TransmissionInequalities(QP.A,Organize);
[QP.Aeq,QP.A,H,QP.f,Organize] = BuildingConstraints(QP.Aeq,QP.A,H,QP.f,Organize,nG,nL,nB);
QP.A = HydroInequalities(QP.A,Organize);
[QP.Aeq,QP.beq,QP.A] = SpinReserveConstraints(QP.Aeq,QP.beq,QP.A,QP.ub,Organize);

%Convert to sparse form
QP.H = spalloc(Organize.t1States,Organize.t1States,Organize.t1States);
for i = 1:1:Organize.t1States
    QP.H(i,i) = H(i);
end
QP.Aeq = sparse(QP.Aeq);
QP.A = sparse(QP.A);
QP.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later
end%Ends function buildMatrices1Step

function [QP,Organize,lb,ub] = countStates1Step(QP,Organize,nG,nB,nL)
%% Organize the order of states 
% states for each generator/storage
% spinning reserve for generator/storage if option is selected
% states for each transmission line/pipe 
% state for P_loss in heating energy balance
% state for spinning reserve shortfall (target - cumulative of all active generators) and SR provided as ancillary service
global Plant
networkNames = fieldnames(Plant.subNet);
xL = 0;lb =[]; ub = [];
Organize.SpinReserveStates = zeros(nG+2,1);
for i = 1:1:nG
    Gen = Plant.Generator(i).QPform;
    if ~isempty(Gen.states)
        [~,fit] = size(Gen.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
        if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage'})
            if ismember(Plant.Generator(i).Type,{'Hydro Storage'})
                s = 1;% Storage treated as generator with 1 state
                lb(end+1,1) = Gen.X.lb; 
                ub(end+1,1) = Gen.X.ub; % Upper bound is MaxGenCapacity
            else
                s = 2;% Storage treated as generator with 2 states (addtl power output relative to 1st dispatch & charging penalty)
                lb(end+1:end+2,1) = [-Gen.Ramp.b(1);0];
                ub(end+1:end+2,1) = [Gen.Ramp.b(2);(1/Gen.Stor.ChargeEff-Gen.Stor.DischEff)*Gen.Ramp.b(2)];
            end 
        else %dispatchable generator or utility
            states = Gen.states(1:nnz(~cellfun('isempty',Gen.states(:,fit))),fit);
            s = length(states);%generator with multiple states
            for k = 1:1:s
                lb(end+1,1) = Gen.(states{k}).lb(fit);
                ub(end+1,1) = Gen.(states{k}).ub(fit);
            end
            if ismember(Plant.Generator(i).Type,{'Electric Generator';'CHP Generator';'Heater';'Chiller';'Absorption Chiller';})
                if isfield(Plant.Generator(i).QPform,'constCost') 
                     QP.constCost(i) = Plant.Generator(i).QPform.constCost;
                end
                if isfield(Plant.Generator(i).VariableStruct,'StartCost') && Plant.Generator(i).VariableStruct.StartCost>0
                    Organize.Dispatchable(i) = 1;
                elseif isfield(Plant.Generator(i).QPform,'constDemand')
                    Organize.Dispatchable(i) = 1;
                elseif QP.constCost(i)>0 || Plant.Generator(i).QPform.(Plant.Generator(i).QPform.states{1,end}).lb(end)>0 % first state has a non-zero lower bound for second optimization or a nonzero cost at an output of zero
                    Organize.Dispatchable(i) = 1;
                end
            end
        end
        QP.organize{i} = xL+1:xL+s; %output is sum of multiple states at each time step
        Organize.States(i)= {xL+1:xL+s};
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') && s ==2
            Organize.Out_vs_State{1,i} = [1,-1];
        elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage'})
            QP.organize{i} = xL+1;  %output state for storage is only SOC (Power in the case of Hydro)
            Organize.Out_vs_State{1,i} = 1;
        else
            Organize.Out_vs_State{1,i} = ones(1,s);
        end
        if Plant.optimoptions.SpinReserve
            include = {'Electric Generator';'CHP Generator';'Hydro Storage';};
            if ismember(Plant.Generator(i).Type,include)
                Organize.SpinReserveStates(i) = xL + s +1; %state of spinning reserve at time 1
                lb(end+1,1) = 0;
                ub(end+1,1) = Plant.Generator(i).Size;
                if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                    lb(end+1,1) = 0;
                    ub(end+1,1) = Plant.Generator(i).VariableStruct.MaxGenCapacity;
                end
            end
            include = {'Electric Storage'};
            if ismember(Plant.Generator(i).Type,include)
                Organize.SpinReserveStates(i) = xL +s+1; %state of spinning reserve at time 1
                lb(end+1,1) = -Plant.Generator(i).QPform.Stor.PeakCharge;
                ub(end+1,1) = Plant.Generator(i).QPform.Stor.PeakDisch;
            end
            xL = xL + s + 1;
        else
            xL = xL + s;
        end
    end
end

%states for transmission lines
for net = 1:1:length(networkNames)
    if ~isempty(Plant.subNet.(networkNames{net}).lineNames)
        limit = Plant.subNet.(networkNames{net}).lineLimit;
        if strcmp(networkNames{net},'Hydro')
            eff = [];
            minimum = Plant.subNet.(networkNames{net}).lineMinimum; 
        else
            eff = Plant.subNet.(networkNames{net}).lineEff;
            minimum = zeros(length(limit(:,1)),1);
        end
        for i = 1:1:length(limit(:,1))
            line  = Plant.subNet.(networkNames{net}).lineNumber(i);
            QP.organize{nG+line} = xL+1; %line state organized into matrix of time vs. line state
            Organize.Out_vs_State{1,nG+line} = 1;
            if isempty(eff) || length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, or river segment, 1 state for each line
                Organize.States(nG+line)= {xL+1};
                lb(end+1,1) = minimum(i);
                ub(end+1,1) = limit(i,1);
                xL = xL + 1;
            else% bi-directional transfer, 3 states for each line (state of the line and penalty term in each direction)
                Organize.States(nG+line)= {[xL+1, xL+2, xL+3]};
                lb(end+1:end+3,1) = [-limit(i,2),0,0];
                ub(end+1:end+3,1) = [limit(i,1), limit(i,1)*(1-eff(i,1)), limit(i,2)*(1-eff(i,2))];
                xL = xL + 3;
            end
        end
    end
end

if Plant.optimoptions.SpinReserve && any(Organize.SpinReserveStates(1:nG)>0)
    Organize.SpinReserveStates(nG+1) = xL+1; %spinning reserve shortfall at t=1
    Organize.SpinReserveStates(nG+2) = xL+2; %SR provided as ancillary service at t=1
    xL = xL + 2; % Add single state for reserve shortfall state which is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    lb(end+1:end+2,1) = 0;
    ub(end+1:end+2,1) = sum(ub(nonzeros(Organize.SpinReserveStates(1,1:nG))));
end

QP.excessHeat = Plant.optimoptions.excessHeat;
QP.excessCool = Plant.optimoptions.excessCool;
[Organize.HeatVented,xL,lb,ub] = VentedEnergy(xL,lb,ub,'H');
[Organize.CoolVented,xL,lb,ub] = VentedEnergy(xL,lb,ub,'C');

for i = 1:1:nB %all buildings have an initial temperature state
    Organize.States(nG+nL+i) = {[xL+1, xL+2, xL+3, xL+4, xL+5]};
    QP.organize{nG+nL+i} = xL+1; %output state for building is temperature
    Organize.Out_vs_State{1,nG+nL+i} = 1;
    xL = xL+5; % 5 states are: Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
    lb(end+1:end+5,1) = [10; 0; 0; 0; 0;];%Dont let building below 10 degrees C
    ub(end+1:end+5,1) = [35; 1e5; 1e5; 10; 10;];%Building stays below 35 degrees C, and excess temperature is less than 10 degrees
end
Organize.t1States = xL;
end%Ends function countStates1Step

function Organize = equaltiyConstraints1Step(Organize)
%% organize equality equations 
% Electric energy balance @ each Electric subNet node  at t = 1
% Heat balance @ each DistricHeat subNet node at t = 1
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
global Plant
networkNames = fieldnames(Plant.subNet);
nG = length(Plant.Generator);
req = 0; % row index of the Aeq matrix and beq vector

% The following puts together the energy balance equations
% 1 energy/mass balance equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Hydro')
        %don't do a water balance, since it depends on multiple time steps.
        %Any extra power at this time step is subtracted from the expected
        %power at the next step (same SOC and flows up river).
    else
        n = length(Plant.subNet.(networkNames{net}).nodes);
        Organize.Balance.(networkNames{net}) = zeros(n,1);
        Organize.Demand.(networkNames{net}) = cell(n,1);
        for i = 1:1:n
            req = req+1;%there is an energy/mass balance at this node
            Organize.Balance.(networkNames{net})(i) = req;
            %%note any demands at this node
            if ~isempty(Plant.subNet.(networkNames{net}).Load{i})
                Organize.Demand.(networkNames{net})(i) = Plant.subNet.(networkNames{net}).Load(i);
            end
        end
    end
end

Organize.Equalities = zeros(nG,2);
Organize.SpinRow = zeros(nG,1);
include = {'Electric Generator';'CHP Generator';'Electric Storage';'Hydro Storage';};
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        %dont link the Power, spill, and outflow with an equality, will
        %link power and outflow with inequality instead
    else
        %link is a field if there is more than one state and the states are linked by an inequality or an equality
        if isfield(Plant.Generator(i).QPform,'link') && isfield(Plant.Generator(i).QPform.link,'eq')
            [m,~] = size(Plant.Generator(i).QPform.link.eq);
            Organize.Equalities(i,1:2) = [req+1, m]; 
            req = req+m;
        end
    end
    if ismember('Electrical',networkNames) && Plant.optimoptions.SpinReserve && ismember(Plant.Generator(i).Type,include)%spinning reserve: no ramping constraint, so SR can be an equalities 
        Organize.SpinRow(i) = req+1;%2 constraints i) ramping, 2) capacity
        req = req+1;
    end
end
Organize.t1Balances = req;
end%Ends function equalityConstraints1Step

function Organize = inequalityConstraints1Step(Organize,nG,nB)
%% Organize inequalities
% No ramping constraints or energy storage inequalities (energy storage is a generator)
% 2 constraints for each bi-directional transmission line
global Plant
networkNames = fieldnames(Plant.subNet);
r = 0; % row index of the A matrix & b vector

%charging penalty eqn
Organize.Inequalities = zeros(nG,1);
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        Organize.Inequalities(i) = r+1;
        r = r+1;
    end
end
%cumulative spinning reserve shortfall: currently only implemented for electric power
if ismember('Electrical',networkNames) && Plant.optimoptions.SpinReserve
    Organize.SpinReserve = r+1;
    r = r+1;
end

%Inequality for Hydro
if ismember('Hydro',networkNames)
    Organize.HydroInequalities = zeros(nG,1);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'Hydro Storage'})
            Organize.HydroInequalities(i) = r+1;
            r = r+1;
        end
    end
end

%Transmission line inequalities (penalty terms)
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        if isfield(Plant.subNet.(networkNames{net}),'lineEff') && ~isempty(Plant.subNet.(networkNames{net}).lineEff)
            eff = Plant.subNet.(networkNames{net}).lineEff;
            for i = 1:1:length(eff(:,1))
                if length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
                    %do nothing, no inequalities linking penalty states
                else%bi-directional power transfer
                    Organize.Transmission(Plant.subNet.(networkNames{net}).lineNumber(i)) = r+1; 
                    r = r+2;
                end
            end
        end
    end
end

for i = 1:1:nB
    Organize.Building.r(i) = r+1;
    r = r+4; % (4 * # of zones) inequalities per building
end
Organize.t1ineq = r;
end%Ends function inewualityConstraints1Step

function [Aeq,beq,H,f] = EqualityConstraints(Aeq,beq,H,f,Organize)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
global Plant
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Hydro')
        %no Hydro mass balance
    else
        if strcmp(networkNames{net},'Electrical')
            out = 'E';
        elseif strcmp(networkNames{net},'DistrictHeat')
            out = 'H';
        elseif strcmp(networkNames{net},'DistrictCool')
            out = 'C';
        end
        for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
            req = Organize.Balance.(networkNames{net})(n);
            equip = Plant.subNet.(networkNames{net}).Equipment{n};
            for k = 1:1:length(equip)
                Gen = Plant.Generator(equip(k)).QPform;
                if ismember(Plant.Generator(equip(k)).Type,{'Hydro Storage'}) 
                    if isfield(Plant.Generator(equip(k)).QPform.output,out)
                        states = Organize.States{equip(k)};
                        Aeq(req,states(1)) = Plant.Generator(equip(k)).QPform.output.(out)(1);
                    end
                elseif ismember(Plant.Generator(equip(k)).Type,{'Electric Storage';'Thermal Storage';}) 
                    if isfield(Plant.Generator(equip(k)).QPform.output,out)
                        states = Organize.States{equip(k)};
                        Aeq(req,states(1)) = Plant.Generator(equip(k)).QPform.output.(out)(1);%additional power beyond 1st dispatch
                        Aeq(req,states(1)+1) = -1;%charging penalty
                    end
                else
                    states = Organize.States{equip(k)};
                    %link is a field if there is more than one state and the states are linked by an inequality or an equality
                    if isfield(Gen,'link') && isfield(Gen.link,'eq')
                        req2 = Organize.Equalities(equip(k),1);
                        for j = 1:1:length(req2)
                            Aeq(req2+(j-1),states) = Gen.link.eq(j,:);
                            beq(req2+(j-1)) = Gen.link.beq(j);
                        end
                    end
                    if ~isempty(Gen.states)
                        [~,fit] = size(Gen.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
                        stateNames = Gen.states(1:nnz(~cellfun('isempty',Gen.states(:,fit))),fit);
                        for j = 1:1:length(stateNames)
                            H(states(j)) = Gen.(stateNames{j}).H(fit);
                            f(states(j)) = Gen.(stateNames{j}).f(fit);
                        end
                        if length(Plant.Generator(equip(k)).QPform.output.(out)(1,:))>1
                            output = Plant.Generator(equip(k)).QPform.output.(out)(:,fit);
                        else output = Plant.Generator(equip(k)).QPform.output.(out);
                        end
                        Aeq(req,states) = output;
                    end
                end
            end
            %%any heat loss term to balance equality
            if strcmp('DistrictHeat',networkNames{net}) && Plant.optimoptions.excessHeat == 1 && Organize.HeatVented(n)~=0
                Aeq(req,Organize.HeatVented(n)) = -1;
            end
        end
    end
end
end%Ends function GeneratorEqualities

function Aeq = LineEqualities(Aeq,Organize)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
global Plant
nG = length(Plant.Generator);
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro') %hydro is skipped because of the time of transfer of the lines
        nodeNames = cell(length(Plant.subNet.(networkNames{net}).nodes),1);
        for n = 1:1:length(nodeNames)
            agregatedNodes = Plant.subNet.(networkNames{net}).nodes{n};
            nodeNames(n) = agregatedNodes(1);
        end
        for k = 1:1:length(Plant.subNet.(networkNames{net}).lineNames)
            line = Plant.subNet.(networkNames{net}).lineNames{k};
            linenum = Plant.subNet.(networkNames{net}).lineNumber(k);
            r = strfind(line,'_');
            origin = line(1:r(1)-1);
            destination = line(r(2)+1:end);
            n1 = find(strcmp(origin,nodeNames),1,'first');
            linestates = Organize.States{nG+linenum};
            req = Organize.Balance.(networkNames{net})(n1);
            if length(linestates)==1 %uni-directional transfer
                Aeq(req,linestates) = -1;
            else
                Aeq(req,linestates) = [-1,0,-1]; %sending node (this node, i)-->(connected node), is positive, thus positive transmission is power leaving the node, the penalty from b->a is power not seen at a
            end
            n2 = find(strcmp(destination,nodeNames),1,'first');
            req = Organize.Balance.(networkNames{net})(n2);
            if length(linestates)==1 %uni-directional transfer
                Aeq(req,linestates) = Plant.subNet.(networkNames{net}).lineEff(k);
            else
                Aeq(req,linestates) = [1,-1,0];%receiving node (connected node)-->(this node, i), is positive, thus positive power is power entering the node, the penalty from a->b is power not seen at b
            end
        end
    end
end
end%Ends function LineEqualities

function A = StorageInequalities(A,Organize)
global Plant;
nG = length(Plant.Generator);
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        states = Organize.States{i};
        r = Organize.Inequalities(i);
        A(r,states(1)) = -(1/Plant.Generator(i).QPform.Stor.ChargeEff-Plant.Generator(i).QPform.Stor.DischEff);
        A(r,states(2)) = -1;
    end
end
end%ends function StorageInequalities

function [Aeq,beq,A] = SpinReserveConstraints(Aeq,beq,A,ub,Organize)
%spinning reserve equalities: individual spinning reserves
%currently only implemented for electric power
global Plant;
if ismember('Electrical',fieldnames(Plant.subNet)) && Plant.optimoptions.SpinReserve
    nG = length(Plant.Generator);
    %individual spinning reserve at each time step (limited by ramp rate and peak gen capacity)
    include = {'Electric Generator';'CHP Generator';'Hydro Storage';'Electric Storage';};
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,include)
            states = Organize.States{i};
            if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                states = states(1);%power output only associated with 1st state
            end
            req = Organize.SpinRow(i);
            Aeq(req+1,states) = 1;
            beq(req+1) = sum(ub(states(1:end-1)));%max power
        end
    end
    A(Organize.SpinReserve,nonzeros(Organize.SpinReserveStates(1,1:nG+1))) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
    A(Organize.SpinReserve,nonzeros(Organize.SpinReserveStates(1,nG+2))) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
end
end%Ends function SpinReserveConstraints

function A = TransmissionInequalities(A,Organize)
%Transmission
%%no connection to previous or later time steps, and no dependence on step size. 
global Plant
nG = length(Plant.Generator);
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if isfield(Plant.subNet.(networkNames{net}),'lineEff')&& ~isempty(Plant.subNet.(networkNames{net}).lineEff)
        eff = Plant.subNet.(networkNames{net}).lineEff;
        for i = 1:1:length(eff(:,1))
            line = Plant.subNet.(networkNames{net}).lineNumber(i);
            lineRow = Organize.Transmission(line);
            states = Organize.States{nG+line};
            if~isempty(lineRow)
                A(lineRow,states) = [(1-eff(i,1)), -1, 0];% Pab*(1-efficiency) < penalty a to b
                A(lineRow+1,states) = [-(1-eff(i,2)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
            end
        end
    end
end
end%Ends function TransmissionInequalities

function A = HydroInequalities(A,Organize)
global Plant
if isfield(Plant.subNet,'Hydro')
    nG = length(Plant.Generator);
    upriver = {};
    upLines = [];
    for i = 1:1:length(Plant.subNet.Hydro.lineNames)
        name = Plant.subNet.Hydro.lineNames{i};
        k = strfind(name,'_');
        upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
        upLines(end+1) = i;
    end
    %% inequality:  Excess Outflow > Power/(Hd*eff*84.67) - Outflow(from 1st dispatch)
    for n = 1:1:length(Plant.subNet.Hydro.nodes)
        %run through the nodes in the hydro network. Generally each node would
        %have only 1 reservior
        equip = Plant.subNet.Hydro.Equipment{n};
        J = upLines(strcmp(Plant.subNet.Hydro.nodes{n},upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
        lineOut = Plant.subNet.Hydro.lineNumber(J);
        for k = 1:1:length(equip)
            if strcmp(Plant.Generator(equip(k)).Type,'Hydro Storage')
                r = Organize.HydroInequalities(equip(k));
                states = Organize.States{equip(k)};
                A(r,states) = Plant.Generator(equip(k)).QPform.Stor.Power2Flow;
                A(r,Organize.States{nG+lineOut})= -1; %Qdownriver (in excess of original solution from multi-time step),  
            end
        end
    end
end
end%Ends function HydroInequalities

function [Aeq,A,H,f,Organize] = BuildingConstraints(Aeq,A,H,f,Organize,nG,nL,nB)
global Plant
for i = 1:1:nB
    states = Organize.States{nG+nL+i};
    %Electric equality %put demands into electrical (heating and cooling) balances
    req = Organize.Balance.Electrical(Plant.Building(i).QPform.nodeE);
    Organize.Building.Electrical.req(i) = req;
    Aeq(req,states(2)) = -Plant.Building(i).QPform.H2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
    Aeq(req,states(3)) = -Plant.Building(i).QPform.C2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
    %Heating equality
%     if Plant.Building(i).QPform.Heating %put into heating eqn
        req = Organize.Balance.DistrictHeat(Plant.Building(i).QPform.nodeH);
        Organize.Building.DistrictHeat.req(i) = req;
        Aeq(req,states(2)) = -1;%subtract building heating needs from heating energy balance
%     end
    %Cooling equality
%     if Plant.Building(i).QPform.Cooling%put into cooling eqn
        req = Organize.Balance.DistrictCool(Plant.Building(i).QPform.nodeC);
        Organize.Building.DistrictCool.req(i) = req;
        Aeq(req,states(3)) = -1;%subtract building cooling needs from cooling energy balance
%     end
    r = Organize.Building.r(i);
    %heating inequality % H = H0 + Hbar where H_bar>=  UA*(Ti-Tset_H) + Cap*(Ti - T(i-1))/dt where dt is in seconds
    % Done in update1Step because of dt: % A(r,states(1)) = (Plant.Building(i).QPform.UA+Plant.Building(i).QPform.Cap*dt);
    A(r,states(2)) = -1;
    %Cooling inequality % Cooling = C0 + C_bar where C_bar>= UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
    % Done in update1Step because of dt: % A(r+1,states(1)) = (Plant.Building(i).QPform.UA-Plant.Building(i).QPform.Cap*dt);
    A(r+1,states(3)) = -1;
    %Upper bound inequality
    A(r+2,states(1)) = 1;% Upper buffer >= Ti - (Tset + Comfort width/2)
    A(r+2,states(4)) = -1;% Upper buffer >= Ti - (Tset + Comfort width/2)
    %Lower bound inequality
    A(r+3,states(1)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
    A(r+3,states(5)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
    %Cost penalty for exceeding temperature bounds
    H(states(4)) = Plant.Building(i).QPform.Discomfort;
    H(states(5)) = Plant.Building(i).QPform.Discomfort;
    f(states(4)) = Plant.Building(i).QPform.Discomfort;
    f(states(5)) = Plant.Building(i).QPform.Discomfort;
end
end%Ends function BuildingConstraints