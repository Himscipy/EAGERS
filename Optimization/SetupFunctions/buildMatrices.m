function QP = buildMatrices(Op,dt)
%builds constant matrices for multi-time-step optimization
%Demands, initial conditions, and utility costs must updated prior to optimization
global Plant 
nG = length(Plant.Generator);
nS = length(dt);
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

Organize.IC = zeros(nG+nL+nB,1);
Organize.States=cell(1,nG+nL+nB);
Organize.Dispatchable = zeros(1,nG);
Organize.Transmission = zeros(nL,1);
Organize.Hydro = false(1,nG);
Organize.Building.r = zeros(nB,1);
Organize.Building.req = zeros(nB,1);
Organize.dt = dt;
Organize.Out_vs_State = cell(1,nG+nL+nB);
QP.organize = cell(nS+1,nG+nL+nB);
QP.constCost = zeros(1,nG);

[QP,Organize,ic] = countIC(QP,Organize,nG,nB,nL);
[QP,Organize,lb,ub] = countStates(QP,Organize,ic,Op,nG,nB,nL);    
totalStates = nS*Organize.t1States+ic;
for t= 2:nS
    for n = 1:1:length(QP.organize(2,:))
        QP.organize{t+1,n} = QP.organize{t,n} + Organize.t1States;
    end
end
H = zeros(totalStates,1);
QP.f = zeros(totalStates,1);
QP.lb = zeros(totalStates,1);
QP.ub = zeros(totalStates,1);

for t= 1:nS
    QP.lb(ic+(t-1)*Organize.t1States+1:ic+t*Organize.t1States)  = lb;
    QP.ub(ic+(t-1)*Organize.t1States+1:ic+t*Organize.t1States)  = ub;
end

[Organize,ec] = equaltiyConstraints(Organize,ic,nG);
totalEqualities = nS*Organize.t1Balances + ic + length(ec);
QP.Aeq = spalloc(totalEqualities,totalStates,totalEqualities*Organize.t1States);
QP.beq = zeros(totalEqualities,1);

Organize = inequalityConstraints(Organize,nG,nB);
totalInequalities = nS*Organize.t1ineq;
QP.A = spalloc(totalInequalities,totalStates,totalInequalities*Organize.t1States);
QP.b = zeros(totalInequalities,1);

%% call sub functions to fill in these matrices now that they are sized
QP.Aeq(1:ic,1:ic) = eye(ic); %initial condition identity matrix
[QP.Aeq,QP.beq,H,QP.f] = GeneratorEqualities(QP.Aeq,QP.beq,H,QP.f,Organize,dt,Op);
QP.Aeq = LineEqualities(QP.Aeq,Organize,dt);
[QP.Aeq,QP.beq] = endStateConstraint(QP.Aeq,QP.beq,QP.organize,ec);
[QP.A,QP.b] = GeneratorInequalities(QP.A,QP.b,Organize,dt);
QP.A = TransmissionInequalities(QP.A,Organize,dt);
[QP.A,QP.b] = SpinReserveInequalities(QP.A,QP.b,Organize,dt);
[QP.Aeq,QP.A,H,QP.f,Organize] = BuildingConstraints(QP.Aeq,QP.A,H,QP.f,Organize,dt,nG,nB,nL,nS);

%convert to sparse matrices
QP.H = spalloc(totalStates,totalStates,totalStates);
for i = 1:1:totalStates
    QP.H(i,i) = H(i);
end
QP.Aeq = sparse(QP.Aeq);
QP.A = sparse(QP.A);
QP.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later
end%Ends function build matrices

function [QP,Organize,ic] = countIC(QP,Organize,nG,nB,nL)
global Plant
% IC for each generator & storage
ic = 0; % row index of the initial condition constraint in the Aeq matrix and beq vector
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Ramp') 
        ic = ic+1;%initial condition state
        QP.organize{1,i} = ic; %output state organized into matrix of time vs. generator (IC)   
        Organize.IC(i) = ic;
    end
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        Organize.Hydro(i) = true;
        ic = ic+1; %second IC for the SOC of the reservior (1st was for Power output)
    end
end
for i = 1:1:nB %all buildings have an initial temperature state
    ic = ic+1;%initial condition state
    QP.organize{1,nG+nL+i} = ic; %output state organized into matrix of time vs. generator (IC)   
    Organize.IC(nG+nL+i) = ic;
end
end%Ends function countIC

function [QP,Organize,lb,ub] = countStates(QP,Organize,ic,Op,nG,nB,nL)
%% organize the order of non-initial condition states 
% states for each generator/storage at t = 1
% spinning reserve for generator/storage at t = 1 if option is selected
% states for each transmission line/pipe at t = 1
% state for P_loss in heating & cooling energy balances
% state for spinning reserve shortfall (target - cumulative of all active generators) and SR provided as ancillary service
% repeat order of generators and lines for t = 2:nS
global Plant
xL = ic; lb =[]; ub = [];
Organize.SpinReserveStates = zeros(nG+2,1);
Organize.Fit = Op;
for i = 1:1:nG
    Gen = Plant.Generator(i).QPform;
    if strcmp(Op,'B') 
        [~,fit] = size(Gen.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
    else
        fit = 1;
    end
    if ~isempty(Gen.states)
        states = Gen.states(1:nnz(~cellfun('isempty',Gen.states(:,fit))),fit);
        s = length(states);%generator with multiple states
        QP.organize{2,i} = xL+1:xL+s; %output is sum of multiple states at each time step
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') && s ==2
            Organize.Out_vs_State{1,i} = [1,-1];
        elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
            QP.organize{2,i} = xL+1; %output state for storage is only SOC (Power in the case of Hydro)
            Organize.Out_vs_State{1,i} = 1;
        else
            Organize.Out_vs_State{1,i} = ones(1,s);
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
        for k = 1:1:s
            lb(end+1,1) = Gen.(states{k}).lb(fit);
            ub(end+1,1) = Gen.(states{k}).ub(fit);
        end
        Organize.States(i)= {xL+1:xL+s};
        xL = xL + s;
        [lb_s,ub_s] = spinReserveGen(i);
        if ~isempty(lb_s)
            Organize.SpinReserveStates(i) = xL+1; %state of spinning reserve at time 1
            lb(end+1,1) = lb_s;
            ub(end+1,1) = ub_s;
            xL = xL + 1;
        end
    end
end

[QP,Organize,xL,lb,ub] = countTransLines(QP,Organize,xL,lb,ub);
if Plant.optimoptions.SpinReserve && any(Organize.SpinReserveStates(1:nG)>0)
    Organize.SpinReserveStates(nG+1) = xL+1; %spinning reserve shortfall at t=1
    Organize.SpinReserveStates(nG+2) = xL+2; %SR provided as ancillary service at t=1
    xL = xL + 2; % Add states for reserve shortfall and ancillary service: shortfall state is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    lb(end+1:end+2,1) = 0;
    ub(end+1:end+2,1) = sum(ub(nonzeros(Organize.SpinReserveStates(1:nG))));
end
[Organize.HeatVented,xL,lb,ub] = VentedEnergy(xL,lb,ub,'H');
[Organize.CoolVented,xL,lb,ub] = VentedEnergy(xL,lb,ub,'C');
for i = 1:1:nB %all buildings have an initial temperature state
    Organize.States(nG+nL+i) = {[xL+1, xL+2, xL+3, xL+4, xL+5]};
    QP.organize{2,nG+nL+i} = xL+1; %output state for building is temperature
    Organize.Out_vs_State{1,nG+nL+i} = 1;
    xL = xL+5; % 5 states are: Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
    lb(end+1:end+5,1) = [10; 0; 0; 0; 0;];%Dont let building below 10 degrees C
    ub(end+1:end+5,1) = [35; 1e5; 1e5; 10; 10;];%Building stays below 35 degrees C, and excess temperature is less than 10 degrees
end
%% Expand by the number of time steps
% # of states per time step is xL 
% Order of states will repeat for nS time steps
Organize.t1States = xL-ic;
end%Ends function countStates

function [lb,ub] = spinReserveGen(i)
global Plant
lb = [];
ub = [];
if Plant.optimoptions.SpinReserve
    include = {'Electric Generator';'CHP Generator';'Hydro Storage';'Electric Storage';};
    if ismember(Plant.Generator(i).Type,include)
        if strcmp(Plant.Generator(i).Type,'Hydro Storage')
            lb(end+1,1) = 0;
            ub(end+1,1) = Plant.Generator(i).VariableStruct.MaxGenCapacity;
        elseif strcmp(Plant.Generator(i).Type,'Electric Storage')
            lb(end+1,1) = -Plant.Generator(i).QPform.Stor.PeakCharge;
            ub(end+1,1) = Plant.Generator(i).QPform.Stor.PeakDisch;
        else
            lb(end+1,1) = 0;
            ub(end+1,1) = Plant.Generator(i).Size;
        end
    end
end
end%Ends function spinReserve

function [QP,Organize,xL,lb,ub] = countTransLines(QP,Organize,xL,lb,ub)
%states for transmission lines
global Plant
networkNames = fieldnames(Plant.subNet);
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
            QP.organize{2,nG+line} = xL+1; %line state organized into matrix of time vs. line state
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
end%Ends function countTransLines

function [Organize,ec] = equaltiyConstraints(Organize,ic,nG)
global Plant
%% Equality Constraints
% IC for each generator and storage device
% Electric energy balance @ each Electric subNet node  at t = 1, including transmission lines
% Heat balance @ each DistricHeat subNet node at t = 1...
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
% Repeat power balance equalities and link equalities at t = 2:nS
networkNames = fieldnames(Plant.subNet);
req = ic; % row index of the Aeq matrix and beq vector

% The following puts together the energy balance equations
% 1 energy/mass balance equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
for net = 1:1:length(networkNames)
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

Organize.Equalities = zeros(nG,2);
for i = 1:1:nG
    %link is a field if there is more than one state and the states are linked by an inequality or an equality
    if isfield(Plant.Generator(i).QPform,'link') && isfield(Plant.Generator(i).QPform.link,'eq')
        [m,~] = size(Plant.Generator(i).QPform.link.eq);
        Organize.Equalities(i,1:2) = [req+1, m]; 
        req = req+m;
    end
end

%%Special case of final state-of-charge equal to initial state or fixed value
ec = [];
if isfield(Plant.optimoptions,'endSOC') && ~strcmp(Plant.optimoptions.endSOC,'Flexible')
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'Electric Storage') || strcmp(Plant.Generator(i).Type,'Thermal Storage') || strcmp(Plant.Generator(i).Type,'Hydro Storage')
            ec(end+1) = i;
        end
    end
end

Organize.t1Balances = req - ic;
end%Ends function equalityConstraints

function Organize = inequalityConstraints(Organize,nG,nB)
global Plant
%% Inequality Constraints
% 2 ramping constraints for each generator at each time step
% constraints for each Energy storage 
% 2 constraints for each bi-directional transmission line
% repeat all of the above for each time 2:nS

r = 0; % row index of the A matrix & b vector
networkNames = fieldnames(Plant.subNet);
%agregate spinning reserve shortfall
%currently only implemented for electric power
if ismember('Electrical',fieldnames(Plant.subNet)) && Plant.optimoptions.SpinReserve
    Organize.SpinReserve = r+1;
    r = r+1; %inequality for net spinning reserve (total of all idividual generator spinning reserves
    include = {'Electric Generator';'CHP Generator';'Electric Storage';'Hydro Storage';};
end
    
%Ramping & Generator Inequalities
Organize.Ramping = zeros(nG,1);
Organize.Inequalities = zeros(nG,2);
Organize.SpinRow = zeros(nG,1);
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Ramp') 
        Organize.Ramping(i) = r+1; 
        r = r+2;
    end
    if isfield(Plant.Generator(i).QPform,'link') && isfield(Plant.Generator(i).QPform.link,'ineq')%link is a field if there is more than one state and the states are linked by an inequality or an equality
        [m,~] = size(Plant.Generator(i).QPform.link.ineq);
        Organize.Inequalities(i,1:2) = [r+1, m];
        r = r+m;
    end
    if ismember('Electrical',networkNames) && Plant.optimoptions.SpinReserve && ismember(Plant.Generator(i).Type,include)%spinning reserve inequalities 
        Organize.SpinRow(i) = r+1;%2 constraints i) ramping, 2) capacity
        r = r+2;
    end
end

%Transmission line inequalities (penalty terms)
for net = 1:1:length(networkNames)
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

for i = 1:1:nB
    Organize.Building.r(i) = r+1;
    r = r+4; % (4 * # of zones) inequalities per building
end

% number of generator inequality constraints & energy imbalances at each time step will be r
% there are 2 ramping constraints on each generator/storage
Organize.t1ineq = r;
end%Ends function inequalityConstraints

function [Aeq,beq,H,f] = GeneratorEqualities(Aeq,beq,H,f,Organize,dt,Op)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
global Plant
nS = length(dt);
t1Balances = Organize.t1Balances;
t1States = Organize.t1States;
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
    end
    if strcmp(networkNames{net},'Hydro')
        Aeq = HydroEqualities(Aeq,Organize,dt);
    else
        T = (1:1:nS)';
        for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
            req = Organize.Balance.(networkNames{net})(n);
            equip = Plant.subNet.(networkNames{net}).Equipment{n}; %only things in this list should have an output or input in this network (electric, heating, cooling
            for k = 1:1:length(equip)
                Gen = Plant.Generator(equip(k)).QPform;
                states = Organize.States{equip(k)};
                %link is a field if there is more than one state and the states are linked by an inequality or an equality
                if isfield(Gen,'link') && isfield(Gen.link,'eq')
                    req2 = Organize.Equalities(equip(k),1);
                    for j = 1:1:length(req2)
                        for t = 1:1:nS
                            Aeq((t-1)*t1Balances+req2+(j-1),(t-1)*t1States+states) = Gen.link.eq(j,:);
                        end
                        beq((T-1)*t1Balances+req2+(j-1)) = Gen.link.beq(j);
                    end
                end
                if ismember(Plant.Generator(equip(k)).Type,{'Electric Storage';'Thermal Storage';})
                    eff = Gen.Stor.DischEff;
                    for t = 1:1:nS
                        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(1)) = -eff/dt(t); %SOC at t
                        if ismember('Y',Gen.states)
                            Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(2)) = -1/dt(t); %charging penalty at t
                        end
                        if t==1
                            Aeq(req,Organize.IC(equip(k))) = eff/dt(t);%SOC at IC
                        else
                            Aeq((t-1)*t1Balances+req,(t-2)*t1States+states(1)) = eff/dt(t); %SOC at t-1
                        end
                    end
                elseif ~isempty(states)
                    if strcmp(Op,'B') 
                        [~,fit] = size(Gen.states);% fit = 2 for generators with 2 different piecewise quadratics when Op = 'B'
                    else
                        fit = 1;
                    end
                    stateNames = Gen.states(1:nnz(~cellfun('isempty',Gen.states(:,fit))),fit);
                    for j = 1:1:length(stateNames)
                        H((T-1)*t1States+states(j)) = Gen.(stateNames{j}).H(fit);%generator  cost fit
                        f((T-1)*t1States+states(j)) = Gen.(stateNames{j}).f(fit);
                    end
                    if length(Gen.output.(out)(1,:))>1
                        output = Gen.output.(out)(1:length(stateNames),fit);
                    else
                        output = Gen.output.(out);
                    end
                    for t = 1:1:nS
                        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states) = output';
                    end
                end
            end
            %%any heat loss term to balance equality
            if strcmp('DistrictHeat',networkNames{net}) && Organize.HeatVented(n)~=0
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+Organize.HeatVented(n)) = -1;
                end
            end
            %%any cool loss term to balance equality
            if strcmp('DistrictCool',networkNames{net}) && Organize.CoolVented(n)~=0
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+Organize.CoolVented(n)) = -1;
                end
            end
        end
    end
end
end%Ends function GeneratorEqualities

function Aeq = HydroEqualities(Aeq,Organize,dt)
%This function loads the mass and energy balances of the hydro network
global Plant
nS = length(dt);
nG = length(Plant.Generator);
t1Balances = Organize.t1Balances;
t1States = Organize.t1States;
upriver = {};
downriver = {};
nodeIndex = (1:1:length(Plant.subNet.Hydro.lineNames));
for i = 1:1:length(Plant.subNet.Hydro.lineNames)
    name = Plant.subNet.Hydro.lineNames{i};
    k = strfind(name,'_');
    upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
    downriver(end+1) = {name(k(2)+1:end)};% node names of the downriver node
end
%% Mass balance
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        n = Plant.Generator(i).QPform.Hydro.subnetNode;
        req = Organize.Balance.Hydro(n);
        J = Plant.subNet.Hydro.lineNumber(strcmp(Plant.subNet.Hydro.nodes{n}(1),upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
        Plant.subNet.Hydro.DownRiverSegment(n) = J;
        if any(strcmp(Plant.subNet.Hydro.nodes{n}(1),downriver))
            Plant.subNet.Hydro.UpRiverNodes(n) = {nodeIndex(strcmp(Plant.subNet.Hydro.nodes{n}(1),downriver))};%upstream nodes
            K = Plant.subNet.Hydro.lineNumber(Plant.subNet.Hydro.UpRiverNodes{n});
        else K = []; %no upstream nodes
        end
        req2 = Organize.Equalities(i,1); %equality for converting power to flow (adding spill flow) = outflow
        for t = 1:1:nS
            %river segments flowing into this node
            for j = 1:1:length(K)
                T = Plant.subNet.Hydro.lineTime(Plant.subNet.Hydro.UpRiverNodes{n}(j));
                tt = sum(dt(1:t));
                if tt<T                            
                    %Do nothing; the inflow rate will be updated in updateMatrices
                elseif tt>=T && tt<=T+dt(1)%between initial condition & first step
                    frac = (tt-T)/dt(1);%portion of flow from step 1, remainder from  SourceSink + (1-frac)*Inflow(t=1) : subtracted from beq in update matrices
                    Aeq((t-1)*t1Balances+req,Organize.States{nG+K(j)})= frac; % Qupriver at step 1, 
                    Plant.subNet.Hydro.frac(Plant.subNet.Hydro.UpRiverNodes{n}(j)) = frac;
                else
                    step = 2;
                    while tt>(T+sum(dt(1:step)))
                        step = step+1;
                    end
                    frac = abs((tt-(T+sum(dt(1:step))))/dt(step));%portion of flow from step, remainder from previous step
                    Aeq((t-1)*t1Balances+req,(step-1)*t1States+Organize.States{nG+K(j)})= frac; % Qupriver at t - floor(T)
                    Aeq((t-1)*t1Balances+req,(step-2)*t1States+Organize.States{nG+K(j)})= (1-frac); % Qupriver at t - ceil(T)
                end 
            end 
            %water flow out of the node
            Aeq((t-1)*t1Balances+req,(t-1)*t1States+Organize.States{nG+J})= -1; %Qdownriver,  
            Plant.Generator(i).QPform.DownRiverSegment = J;

            states = Organize.States{i};
            %SOC of the reservior in 1000 acre ft.
            Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(2)) = -12.1/dt(t); %SOC at t in acre-ft / hours converted to 1000 ft^3/s:  12100 ft^3/s * 3600s = 1000 acre ft
            if t==1
                Aeq(req,Organize.IC(i)+1) = 12.1/dt(t);%SOC at IC is 1 after dam power IC
            else
                Aeq((t-1)*t1Balances+req,(t-2)*t1States+states(2)) = 12.1/dt(t); %SOC at t-1 in acre-ft / hours converted to 1000 ft^3/s:  12100 ft^3/s * 3600s = 1000 acre ft
            end
            %% convert power to outflow (other part of balance is done in Generator Equalities when doing the electric network, this is the - Qdownriver
            Aeq((t-1)*t1Balances+req2,(t-1)*t1States+Organize.States{nG+J})= -1; % Power/(eff*head*84.76) + Spill - Qdownriver = 0
        end 
    else
        %% add water district here
    end
end
end%Ends function Hydro Equalities


function Aeq = LineEqualities(Aeq,Organize,dt)
%this function loads generators/chillers, etc into the apropriate energy balance of Aeq
global Plant
nS = length(dt);
nG = length(Plant.Generator);
t1Balances = Organize.t1Balances;
t1States = Organize.t1States;
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro') %hydro is done later because of the time of transfer of the lines
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
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+linestates) = -1;
                end
            else
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+linestates) = [-1,0,-1]; %sending node (this node, i)-->(connected node), is positive, thus positive transmission is power leaving the node, the penalty from b->a is power not seen at a
                end
            end
            n2 = find(strcmp(destination,nodeNames),1,'first');
            req = Organize.Balance.(networkNames{net})(n2);
            if length(linestates)==1 %uni-directional transfer
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+linestates) = Plant.subNet.(networkNames{net}).lineEff(k);
                end
            else
                for t = 1:1:nS
                    Aeq((t-1)*t1Balances+req,(t-1)*t1States+linestates) = [1,-1,0];%receiving node (connected node)-->(this node, i), is positive, thus positive power is power entering the node, the penalty from a->b is power not seen at b
                end
            end
        end
    end
end
end%Ends function LineEqualities

function [A,b] = GeneratorInequalities(A,b,Organize,dt)
global Plant
nG = length(Plant.Generator);
nS = length(dt);
t1ineq = Organize.t1ineq;
t1States = Organize.t1States;
for i = 1:1:nG
    states = Organize.States{i};
    %Ramping 
    if Organize.Ramping(i)>0
        %%if storage, ramping only affects 1st state
        if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage'})
            rampStates = states(1);
        else
            rampStates = states;
        end
        r = Organize.Ramping(i);
        for t = 1:1:nS
            A((t-1)*t1ineq+r,(t-1)*t1States+rampStates) = 1/dt(t);%ramp up 
            A((t-1)*t1ineq+r+1,(t-1)*t1States+rampStates) = -1/dt(t);%ramp down
            if t ==1 %ramping from initial condition
                A(r,Organize.IC(i)) = -1/dt(t); %ramp up 
                A(r+1,Organize.IC(i)) = 1/dt(t); %ramp down 
            else %condition at previous time step
                A((t-1)*t1ineq+r,(t-2)*t1States+rampStates) = -1/dt(t); %ramp up 
                A((t-1)*t1ineq+r+1,(t-2)*t1States+rampStates) = 1/dt(t); %ramp down 
            end
            b((t-1)*t1ineq+r,1) = Plant.Generator(i).QPform.Ramp.b(1); %ramp up 
            b((t-1)*t1ineq+r+1,1) = Plant.Generator(i).QPform.Ramp.b(2); %ramp down 
        end
    end
    %Inequalities constraints
    if Organize.Inequalities(i,1)>0
        ineqRow = Organize.Inequalities(i,1);
        for t = 1:1:nS
            for k = 1:1:Organize.Inequalities(i,2)
                if k == 1 && ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage'}) && ismember('Y',Plant.Generator(i).QPform.states)
                    A((t-1)*t1ineq+ineqRow,(t-1)*t1States+states(1)) = Plant.Generator(i).QPform.link.ineq(1,1);  % SOC at t  
                    A((t-1)*t1ineq+ineqRow,(t-1)*t1States+states(2)) = Plant.Generator(i).QPform.link.ineq(1,2);  % charging state at t (-1)
                    if t ==1 %SOC change from IC
                        A(ineqRow,Organize.IC(i)) = -Plant.Generator(i).QPform.link.ineq(1,1); % SOC at t-1
                    else
                        A((t-1)*t1ineq+ineqRow,(t-2)*t1States+states(1)) = -Plant.Generator(i).QPform.link.ineq(1,1);  % SOC at t-1
                    end
                    loss = (Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
                    b((t-1)*t1ineq+ineqRow+(k-1)) = Plant.Generator(i).QPform.link.bineq(k) - loss*(1-Plant.Generator(i).QPform.Stor.ChargeEff);
                else
                    A((t-1)*t1ineq+ineqRow+(k-1),(t-1)*t1States+states) = Plant.Generator(i).QPform.link.ineq(k,:); %typically buffer on SOC
                    b((t-1)*t1ineq+ineqRow+(k-1)) = Plant.Generator(i).QPform.link.bineq(k);
                end
            end
        end
    end
end
end%Ends function GeneratorInequalities

function A = TransmissionInequalities(A,Organize,dt)
%Transmission
%%no connection to previous or later time steps, and no dependence on step size. 
global Plant
nS = length(dt);
nG = length(Plant.Generator);
t1ineq = Organize.t1ineq;
t1States = Organize.t1States;
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    if isfield(Plant.subNet.(networkNames{net}),'lineEff')&& ~isempty(Plant.subNet.(networkNames{net}).lineEff)
        eff = Plant.subNet.(networkNames{net}).lineEff;
        for i = 1:1:length(eff(:,1))
            line = Plant.subNet.(networkNames{net}).lineNumber(i);
            lineRow = Organize.Transmission(line);
            states = Organize.States{nG+line};
            if~isempty(lineRow)
                for t = 1:1:nS
                    A((t-1)*t1ineq+lineRow,(t-1)*t1States+states) = [(1-eff(i,1)), -1, 0];% Pab*(1-efficiency) < penalty a to b
                    A((t-1)*t1ineq+lineRow+1,(t-1)*t1States+states) = [-(1-eff(i,2)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
                end
            end
        end
    end
end
end%Ends function TransmissionInequalities

function [A,b] = SpinReserveInequalities(A,b,Organize,dt)
%spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
%currently only implemented for electric power
global Plant;
if ismember('Electrical',fieldnames(Plant.subNet)) && Plant.optimoptions.SpinReserve
    nG = length(Plant.Generator);
    nS = length(dt);
    t1ineq = Organize.t1ineq;
    t1States = Organize.t1States;
    %total spinning reserve shortfall at each time step
    SRstates = nonzeros(Organize.SpinReserveStates(1:nG+1));
    SRancillary = Organize.SpinReserveStates(nG+2);
    for t = 1:1:nS
        A((t-1)*t1ineq+Organize.SpinReserve,(t-1)*t1States+SRstates) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        A((t-1)*t1ineq+Organize.SpinReserve,(t-1)*t1States+SRancillary) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
    end
    %individual spinning reserve at each time step (limited by ramp rate and peak gen capacity)
    include = {'Electric Generator';'CHP Generator';'Hydro Storage';};
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,include)
            SRstate = Organize.SpinReserveStates(i);
            states = Organize.States{i};
            if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                states = states(1);%power output only associated with 1st state
            end
            r = Organize.SpinRow(i);
            for t = 1:1:nS
                A((t-1)*t1ineq+r,(t-1)*t1States+SRstate) = 1/dt(t); %SR + power(t) - power(t-1)<= ramprate*dt
                A((t-1)*t1ineq+r,(t-1)*t1States+states) = 1/dt(t);
                if t ==1 %ramping from IC
                    A(r,Organize.IC(i)) = -1/dt(t); % Power at t-1
                else
                    A((t-1)*t1ineq+r,(t-2)*t1States+states) = -1/dt(t); % Power at t-1
                end
                b((t-1)*t1ineq+r) = Plant.Generator(i).QPform.Ramp.b(1);%ramp up constraint

                A((t-1)*t1ineq+r+1,(t-1)*t1States+SRstate) = 1; %SR + power <= Size
                A((t-1)*t1ineq+r+1,(t-1)*t1States+states) = 1;
                b((t-1)*t1ineq+r+1) = Plant.Generator(i).Size; %max capacity constraint
            end
        end
    end
    %individual spinning reserve at each time step (limited by discharge rate and storage capacity)
    include = {'Electric Storage';};
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,include)
            SRstate = Organize.SpinReserveStates(i);
            states = Organize.States{i};
            SOC = states(1);
            r = Organize.SpinRow{i};
            for t = 1:1:nS
                eff = Plant.Generator(i).QPform.Stor.DischEff;
                A((t-1)*t1ineq+r,(t-1)*t1States+SRstate) = 1; %SR + eff*(SOC(t-1) - SOC(t))/dt <= peak discharge
                A((t-1)*t1ineq+r,(t-1)*t1States+SOC) = -eff/dt(t);
                A((t-1)*t1ineq+r+1,(t-1)*t1States+SRstate) = 1; %SR - SOC(t-1)/dt <= 0
                if t ==1 %SOC change from IC
                    A(r,Organize.IC(i)) = eff/dt(t); % SOC at t-1
                    A(r+1,Organize.IC(i)) = -eff/dt(t); % SOC at t-1
                else
                    A((t-1)*t1ineq+r,(t-2)*t1States+SOC) = eff/dt(t); % SOC at t-1
                    A((t-1)*t1ineq+r+1,(t-1)*t1States+SOC) = -eff/dt(t); % SOC at t-1
                end
                b((t-1)*t1ineq+r) = Plant.Generator(i).QPform.Ramp.b(2);%peak discharge constraint
            end
        end
    end
end
end%Ends function SpinReserveInequalities

function [Aeq,beq] = endStateConstraint(Aeq,beq,organize,ec)
global Plant
n = length(beq) - length(ec);
for j = 1:1:length(ec)
    n = n+1;
    i = ec(j);
    if isnumeric(Plant.optimoptions.endSOC)%assume a vector of fixed end value conditions
        Aeq(n,organize{end,i}) = 1;
        beq(n,1) = Plant.optimoptions.endSOC(i);
    elseif strcmp(Plant.optimoptions.endSOC,'Initial')
        beq(n,1) = 0;
        if strcmp(Plant.Generator(i).Type,'Hydro Storage')
            Aeq(n,organize{1,i}+1) = 1;
            Aeq(n,organize{end,i}+1) = -1;
        else%electric and thermal storage
            Aeq(n,organize{1,i}) = 1;
            Aeq(n,organize{end,i}) = -1;
        end
    end
end
end%Ends function endStateConstraint

function [Aeq,A,H,f,Organize] = BuildingConstraints(Aeq,A,H,f,Organize,dt,nG,nB,nL,nS)
global Plant
%States representing building are: Tzone, Hating, Cooling, Excess T, shortfall T
t1Balances = Organize.t1Balances;
t1ineq = Organize.t1ineq;
t1States = Organize.t1States;
for i = 1:1:nB
    states = Organize.States{nG+nL+i};
    %Electric equality
    req = Organize.Balance.Electrical(Plant.Building(i).QPform.nodeE);
    Organize.Building.Electrical.req(i) = req;
    for t = 1:1:nS
        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(2)) = -Plant.Building(i).QPform.H2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(3)) = -Plant.Building(i).QPform.C2E;% Electricity = E0 + H2E*(Heating-H0) + C2E*(Cooling-C0)
    end
    %Heating equality
    req = Organize.Balance.DistrictHeat(Plant.Building(i).QPform.nodeH);
    Organize.Building.DistrictHeat.req(i) = req;
    for t = 1:1:nS
        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(2)) = -1;%subtract building heating needs from heating energy balance
    end
    %Cooling equality
    req = Organize.Balance.DistrictCool(Plant.Building(i).QPform.nodeC);
    Organize.Building.DistrictCool.req(i) = req;
    for t = 1:1:nS
        Aeq((t-1)*t1Balances+req,(t-1)*t1States+states(3)) = -1;%subtract building cooling needs from cooling energy balance
    end
    r = Organize.Building.r(i);

    for t = 1:1:nS
        %heating inequality % H = H0 + Hbar where H_bar>=  UA*(Ti-Tset_H) + Cap*(Ti - T(i-1))/dt where dt is in seconds
        A((t-1)*t1ineq+r,(t-1)*t1States+states(1)) = (Plant.Building(i).QPform.UA+Plant.Building(i).QPform.Cap/(3600*dt(t)));
        A((t-1)*t1ineq+r,(t-1)*t1States+states(2)) = -1;
        if t == 1
            A(r,Organize.IC(nG+nL+i)) = -Plant.Building(i).QPform.Cap/(3600*dt(t));
        else
            A((t-1)*t1ineq+r,(t-2)*t1States+states(1)) = -Plant.Building(i).QPform.Cap/(3600*dt(t));
        end
        %Cooling inequality % Cooling = C0 + C_bar where C_bar>= UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
        A((t-1)*t1ineq+r+1,(t-1)*t1States+states(1)) = -Plant.Building(i).QPform.UA - Plant.Building(i).QPform.Cap/(3600*dt(t));
        A((t-1)*t1ineq+r+1,(t-1)*t1States+states(3)) = -1;
        if t == 1
            A(r+1,Organize.IC(nG+nL+i)) = -Plant.Building(i).QPform.Cap/(3600*dt(t));
        else
            A((t-1)*t1ineq+r+1,(t-2)*t1States+states(1)) = -Plant.Building(i).QPform.Cap/(3600*dt(t));% Cooling>= C0 + UA*(Tset_C-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
        end
        %Upper bound inequality
        A((t-1)*t1ineq+r+2,(t-1)*t1States+states(1)) = 1;% Upper buffer >= Ti - (Tset + Comfort width/2)
        A((t-1)*t1ineq+r+2,(t-1)*t1States+states(4)) = -1;% Upper buffer >= Ti - (Tset + Comfort width/2)
        %Lower bound inequality
        A((t-1)*t1ineq+r+3,(t-1)*t1States+states(1)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
        A((t-1)*t1ineq+r+3,(t-1)*t1States+states(5)) = -1;% Lower buffer >= (Tset - Comfort width/2) - Ti
        %Cost penalty for exceeding temperature bounds
        H((t-1)*t1States+states(4)) = dt(t)*Plant.Building(i).QPform.Discomfort;
        H((t-1)*t1States+states(5)) = dt(t)*Plant.Building(i).QPform.Discomfort;
        f((t-1)*t1States+states(4)) = dt(t)*Plant.Building(i).QPform.Discomfort;
        f((t-1)*t1States+states(5)) = dt(t)*Plant.Building(i).QPform.Discomfort;
    end
end
end%Ends function BuildingConstraints