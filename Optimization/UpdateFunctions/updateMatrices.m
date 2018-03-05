function QP = updateMatrices(QP,Date,scaleCost,marginCost,Forecast,EC)
% QP is the set of optimization matrices 
% IC is the intial condition
% Time is the vector of time steps
% scale cost is the matrix multiplier of cost for each generator at each time step
% MarginCost is the marginal cost of generation (used in the final value of the energy storage)
% Demand is the forecasted loads
% Renewable is the forecasted uncontrollable generation
% EC is the end condition for the threshold optimization
global Plant CurrentState
if strcmp(Plant.optimoptions.solver,'ANN')
    QP.solver = 'quadprog';
else
    QP.solver = Plant.optimoptions.solver;
end
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
networkNames = fieldnames(Plant.subNet);
dt = (Date(2:end) - Date(1:end-1))*24;
nS = length(dt);
nL = length(QP.Organize.IC)-nG-nB;
%% update demands and storage self-discharge
for net = 1:1:length(networkNames)
    out = Plant.subNet.(networkNames{net}).abbreviation;
    if strcmp(networkNames{net},'Electrical')
        if Plant.optimoptions.SpinReserve
            QP.b(QP.Organize.SpinReserve) = -Forecast.SRtarget;% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        end
        if isfield(Forecast,'Renewable')
            QP.Renewable = Forecast.Renewable;
        end
    end
    QP.constDemand.(out).req = zeros(length(Plant.subNet.(networkNames{net}).nodes)*nS,1);
    QP.constDemand.(out).load = zeros(length(Plant.subNet.(networkNames{net}).nodes)*nS,nG);
    for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes) %run through all the nodes in this network
        equip = Plant.subNet.(networkNames{net}).Equipment{n}; %equipment at this node
        req = QP.Organize.Balance.(networkNames{net})(n);%balance at this node (t = 1)
        req = req:QP.Organize.t1Balances:((nS-1)*QP.Organize.t1Balances+req);%balance at this node (t = 1:nS)
        load = nonzeros(QP.Organize.Demand.(networkNames{net}){n}); %loads at this node
        if ~isempty(load) %need this in case there is no field Forecast.Demand
            QP.beq(req) = sum(Forecast.Demand.(out)(:,load),2); %multiple demands can be at the same node, or none
        end
        if isfield(Forecast,'Renewable')% subtract renewable generation 
            for k = 1:1:length(equip)
                i = equip(k);
                if strcmp(Plant.Generator(i).Source,'Renewable')
                    if strcmp(Plant.Generator(i).Type,'Solar') && strcmp(networkNames{net},'DirectCurrent')
                        QP.beq(req) = QP.beq(req) - Forecast.Renewable(:,i); %put renewable generation into energy balance at correct node
                    end
                end
            end
        end
        for j = 1:1:length(equip)
            k = equip(j);
            if ~isempty(strfind(Plant.Generator(k).Type,'Storage'))
                if isfield(Plant.Generator(k).QPform.output,out)
                    loss = (Plant.Generator(k).QPform.Stor.SelfDischarge*Plant.Generator(k).QPform.Stor.UsableSize*Plant.Generator(k).QPform.Stor.DischEff);
                    QP.beq(req) = QP.beq(req) + loss; %account for self-discharge losses
                end
            end 
            if isfield(Plant.Generator(k).QPform,'constDemand')  && isfield(Plant.Generator(k).QPform.constDemand,out) && strcmp(QP.Organize.Fit,'B') %when using fit B record the constant electrical demands of the chillers
                QP.constDemand.(out).req((n-1)*nS+1:n*nS,1) = req;
                QP.constDemand.(out).load((n-1)*nS+1:n*nS,k) = Plant.Generator(k).QPform.constDemand.(out);
            end
        end
    end
    if strcmp(networkNames{net},'Hydro')
        for n = 1:1:length(Plant.subNet.Hydro.nodes) %run through all the nodes in this network
            %% Node inflows (source/sink terms and upstream flow at time t-T ago if t<T)
            %need to be able to forecast sink/source
            req = QP.Organize.Balance.Hydro(n);%mass balance at this node (t = 1)
            req = req:QP.Organize.t1Balances:((nS-1)*QP.Organize.t1Balances+req);%balance at this node (t = 1:nS)
            QP.beq(req)= QP.beq(req) - Forecast.Hydro.InFlow(:,n); %source/sink term
        end
    end
end  

%Building inequalities
for i = 1:1:nB %all Forecasts are vectors of nS steps
    states = QP.Organize.States{nG+nL+i};
    T_state = states(1):QP.Organize.t1States:(nS-1)*QP.Organize.t1States + states(1);

    r = QP.Organize.Building.r(i):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq + QP.Organize.Building.r(i);
    QP.b(r) = Plant.Building(i).QPform.UA*Forecast.Building.Tset_H(:,i);%heating inequality
    QP.b(r+1) = -Plant.Building(i).QPform.UA*Forecast.Building.Tset_C(:,i);%cooling inequality
    QP.b(r+2) = Forecast.Building.Tmax(:,i);%upper buffer inequality
    QP.b(r+3) = -Forecast.Building.Tmin(:,i);%lower buffer inequality
    
    req = QP.Organize.Building.Electrical.req(i):QP.Organize.t1Balances:(nS-1)*QP.Organize.t1Balances + QP.Organize.Building.Electrical.req(i);
    QP.beq(req) = QP.beq(req) + Forecast.Building.E0(:,i);%electric load
    %% update heating eqn
    req = QP.Organize.Building.DistrictHeat.req(i):QP.Organize.t1Balances:(nS-1)*QP.Organize.t1Balances + QP.Organize.Building.DistrictHeat.req(i);
    QP.Organize.Building.H_Offset(:,i) = Plant.Building(i).QPform.UA*(Forecast.Building.Tzone(:,i)-Forecast.Building.Tset_H(:,i)) - Forecast.Building.H0(:,i);
    QP.beq(req) = QP.beq(req) - QP.Organize.Building.H_Offset(:,i);%Heating load
    QP.lb(T_state+1) = QP.lb(T_state+1) + QP.Organize.Building.H_Offset(:,i);%Ensures net building heating is greater than zero
    
    %% update cooling eqn
    req = QP.Organize.Building.DistrictCool.req(i):QP.Organize.t1Balances:(nS-1)*QP.Organize.t1Balances + QP.Organize.Building.DistrictCool.req(i);
    QP.Organize.Building.C_Offset(:,i) = Plant.Building(i).QPform.UA*(Forecast.Building.Tset_C(:,i)-Forecast.Building.Tzone(:,i)) - Forecast.Building.C0(:,i);
    QP.beq(req) = QP.beq(req) - QP.Organize.Building.C_Offset(:,i);%Cooling load
    QP.lb(T_state+2) = QP.lb(T_state+2) + QP.Organize.Building.C_Offset(:,i);%Ensures net building cooling is greater than zero
end

%Updating lower and upper bound for Hydro resevoir SOC
if isfield(Plant.subNet,'Hydro')
    QP.Organize.hydroSOCoffset = zeros(1,length(Plant.subNet.Hydro.nodes));
    for i = 1:1:nG %Make sure bounds are between 0 and the maximum usable size for each generator
        if ismember(Plant.Generator(i).Type,{'Hydro Storage'})
            n = Plant.Generator(i).QPform.Hydro.subnetNode;%dam #
            if isfield(Plant,'WYForecast') && Date(end)<Plant.WYForecast{end}.Timestamp(end)
                QP.Organize.hydroSOCoffset(n) = interp1(Plant.WYForecast{end}.Timestamp,Plant.WYForecast{end}.hydroSOC(:,n),Date(end));% re-center so zero is actual where this target is
                states = QP.Organize.States{i}; %States for Reservoirs
                SOCstates = states(2):QP.Organize.t1States:(nS-1)*QP.Organize.t1States + states(2); %All time steps for resevoir SOC
                range = max(0.1*Plant.Generator(i).QPform.Stor.UsableSize,1+abs(CurrentState.Hydro(i)-QP.Organize.hydroSOCoffset(n)));
                QP.ub(SOCstates) = min(range,Plant.Generator(i).QPform.Stor.UsableSize-QP.Organize.hydroSOCoffset(n));
                QP.lb(SOCstates) = max(-range,-QP.Organize.hydroSOCoffset(n));
            end
            SOCnow = CurrentState.Hydro(i)-QP.Organize.hydroSOCoffset(n);
            QP.beq(QP.Organize.IC(i)+1) = SOCnow;%SOC in reservior (2nd ic for dam)
            QP.ub(QP.Organize.IC(i)+1) = SOCnow+1;%+1 just to help solver find feasible
            QP.lb(QP.Organize.IC(i)+1) = SOCnow-1;%-1 just to help solver find feasible
        end
    end
end

%Adding cost for Line Transfer Penalties to differentiate from spillway flow
if ismember('Hydro',networkNames) && ~isempty(Plant.subNet.Electrical.lineNames)
    for k = 1:1:length(Plant.subNet.Electrical.lineNames)
        line = Plant.subNet.Electrical.lineNumber(k);
        states = QP.Organize.States{nG+line};
        if length(states)>1 %bi-directional transfer with penalty states
            s2 = states(2):QP.Organize.t1States:((nS-1)*QP.Organize.t1States+states(2));%penalty from a to b
            s3 = states(3):QP.Organize.t1States:((nS-1)*QP.Organize.t1States+states(3));% penalty from b to a
            QP.f(s2) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
            QP.f(s3) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
        end
    end 
end 

%update initial conditions
for i = 1:1:nG+nL+nB
    if QP.Organize.IC(i)>0 
        if i<=nG %generators and storage devices
            QP.beq(QP.Organize.IC(i)) = CurrentState.Generators(i);
            QP.ub(QP.Organize.IC(i)) = CurrentState.Generators(i)+1;%+1 just to help solver find feasible
        else %all buildings have an initial air zone temperature state
            QP.beq(QP.Organize.IC(i)) = CurrentState.Buildings(1,i-nG-nL);
            QP.ub(QP.Organize.IC(i)) = CurrentState.Buildings(1,i-nG-nL)+1;%+1 just to help solver find feasible
        end
    end
end


%update costs
H = diag(QP.H);
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i};
        allStates = (0:QP.Organize.t1States:(nS-1)*QP.Organize.t1States)'*ones(1,length(states))+ones(nS,1)*states;
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') && Plant.Generator(i).VariableStruct.SellBackRate>0
            QP.f(allStates(:,1)) = QP.f(allStates(:,1)).*scaleCost(:,i).*dt;
            if Plant.Generator(i).VariableStruct.SellBackRate == -1
                QP.f(allStates(:,2)) = QP.f(allStates(:,2)).*scaleCost(:,i).*dt; %sellback is a fixed percent of purchase costs (percent set wehn building matrices)
            else
                %constant sellback rate was taken care of when building the matrices
            end
        elseif ismember(Plant.Generator(i).Type,{'Utility';'Electric Generator';'CHP Generator';'Chiller';'Heater';'Hydrogen Generator';'Cooling Tower';'Electrolyzer'})%all generators and utilities (without sellback)
            for j = 1:1:length(states)
                H(allStates(:,j)) = H(allStates(:,j)).*scaleCost(:,i).*dt; 
                QP.f(allStates(:,j)) = QP.f(allStates(:,j)).*scaleCost(:,i).*dt; 
            end          
        elseif isfield(Plant.Generator(i).QPform,'Stor') % storage systems
            if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                s_end = (nS-1)*QP.Organize.t1States + states(2); %SOC is second state
            else
                s_end = (nS-1)*QP.Organize.t1States + states(1);
                if isfield(Plant.Generator(i).QPform,'Y')%charging penalty state
                    QP.f(states(2):QP.Organize.t1States:s_end+1) = 1e-6; %small penalty on charging state so that if there are no other generators it keeps the charging constraint as an equality
                end
            end
            
            StorSize = Plant.Generator(i).QPform.Stor.UsableSize;
            if isfield(Plant.Generator(i).QPform,'U') %has buffer
                BuffSize = Plant.Generator(i).QPform.U.ub;
            else
                BuffSize = 0;
            end
            if isempty(EC) %full forecast optimization
                S = fieldnames(Plant.Generator(i).QPform.output);
                Max = 1.05*marginCost.(S{1}).Max;
                Min = 0.75*marginCost.(S{1}).Min;
                a1 = -Max; % fitting C = a1*SOC + 0.5*a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                a2 = (Max - Min)/(StorSize);
                H(s_end) = a2;%quadratic final value term loaded into SOC(t=nS)  %quadprog its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = a1;%linear final value term loaded into SOC(t=nS)
                if BuffSize>0
                    QP.f(allStates(:,end-1)) = Min;%this is the linear buffer term loaded into the lower & upper buffer
                    H(allStates(:,end-1)) = (2*Max-Min)/BuffSize;%this is the quadratic buffer term loaded into the lower & upper buffer
                    QP.f(allStates(:,end)) = Min;%this is the linear buffer term loaded into the lower & upper buffer
                    H(allStates(:,end)) = (2*Max-Min)/BuffSize;%this is the quadratic buffer term loaded into the lower & upper buffer
                end
            else %used in Online loop when there is a target EC determined by dispatch loop
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
%                 if strcmp(Plant.Generator(i).Type,'Electric Storage')
%                     type = 'Electrical';
%                 elseif strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Heat')
%                     type = 'DistrictHeat';
%                 elseif strcmp(Plant.Generator(i).Type,'Thermal Storage')
%                     type = 'DistrictCool';
%                 elseif strcmp(Plant.Generator(i).Type,'Hydro Storage')
%                     type = 'Hydro';
%                 end
%                 rows = QP.Organize.Inequalities{i};
%                 nr = length(states)/nS;
%                 PeakChargePower = Plant.Generator(i).QPform.Ramp.b(1);
%                 dSOC_10perc = .1*PeakChargePower*(Date(end)-Date(1)); %energy (kWh) if charging at 10%
%                 H(s_end) = -2*marginCost.(type).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
%                 QP.f(s_end) = -marginCost.(type).Min;%linear final value term loaded into SOC(t=nS)
%                 nIC = nnz(QP.Organize.IC(1:i)); %order in IC
%                 QP.beq(nIC) = IC(I)-EC(I);
%                 for t = 1:1:nS
%                     Xn = states((t-1)*nt+1);
%                     Rn = rows(t*nr);
%                     QP.lb(Xn) = -EC(i);%change lb so that SOC = 0 coresponds to EC
%                     QP.ub(Xn) = StorSize - EC(i);%change ub so that SOC = 0 coresponds to EC
%                     QP.b(Rn-1) = -BuffSize + EC(i);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
%                     QP.b(Rn) = StorSize-BuffSize - EC(i);%change lb so that SOC = 0 coresponds to EC
%                 end       
            end
        end
    end
end 
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(:,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    allSRshort = (0:QP.Organize.t1States:(nS-1)*QP.Organize.t1States)'+ones(nS,1)*SRshort;
%     SRancillary = QP.Organize.SpinReserveStates(:,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt./Forecast.SRtarget;% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        req = nonzeros(QP.Organize.Balance.Electrical);
        SpinCost = zeros(nS,1);
        for t = 1:1:nS
            SpinCost(t) = 2*dt(t)./(0.05*sum(QP.beq(req)));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
            req = req + QP.Organize.t1Balances;
        end
    end
    H(allSRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(allSRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);