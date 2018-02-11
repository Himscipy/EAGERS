function QP = updateMatrices1Step(QP,Forecast,scaleCost,marginal,StorPower,dt,IC,FirstProfile,limit,t,T)
% update the equalities with the correct demand, and scale fuel and electric costs
% EC is the expected end condition at this time stamp (can be empty)
% StorPower is the expected output/input of any energy storage at this timestep (can be empty)
% MinPower and MaxPower define the range of this generator at this timestep
%T is the building temperatures
global Plant
QP.solver = Plant.optimoptions.solver;
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
nL = length(QP.Organize.States)-nB-nG;

UB = zeros(1,nG);
dX_dt = zeros(1,nG);
for i = 1:1:nG
    if ~isempty(Plant.Generator(i).QPform.states)
        if strcmp(Plant.Generator(i).Type,'Hydro Storage')
            %%do something to set UB
        else
            states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
            for j = 1:1:length(states)
                UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(end);
            end
        end
    end
    dX_dt(i) = Plant.Generator(i).VariableStruct.dX_dt;
end
if strcmp(limit, 'unconstrained')
    MinPower = 0*UB;
    MaxPower = UB;
elseif strcmp(limit, 'constrained')
    MinPower = max(0,IC(1:nG)-dX_dt*dt(t));
    MaxPower = min(UB,IC(1:nG)+dX_dt*dt(t));
else
    MinPower = max(0,IC(1:nG)-dX_dt*sum(dt(1:t)));
    MaxPower = min(UB,IC(1:nG)+dX_dt*sum(dt(1:t)));
end

networkNames = fieldnames(Plant.subNet);
%% update demands and storage self-discharge
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Hydro')
        %don't do a water balance, since it depends on multiple time steps.
        %Any extra outflow at this time step is subtracted from the expected
        %outflow at the next step (same SOC and flows up river).
    else
        if strcmp(networkNames{net},'Electrical')
            out = 'E';
            if Plant.optimoptions.SpinReserve
                QP.b(QP.Organize.SpinReserve) = -Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E(t,:),2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
            end
        elseif strcmp(networkNames{net},'DistrictHeat')
            out = 'H';
        elseif strcmp(networkNames{net},'DistrictCool')
            out = 'C';
        end
        QP.constDemand.(out).req = zeros(length(Plant.subNet.(networkNames{net}).nodes),1);
        QP.constDemand.(out).load = zeros(length(Plant.subNet.(networkNames{net}).nodes),nG);
        for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
            equip = Plant.subNet.(networkNames{net}).Equipment{n};
            req = QP.Organize.Balance.(networkNames{net})(n);
            loads = nonzeros(QP.Organize.Demand.(networkNames{net}){n});
            if ~isempty(loads)%need this to avoid the next line when there is no load
                QP.beq(req) = sum(Forecast.Demand.(out)(t,loads),2); %multiple demands can be at the same node
            end
            for j = 1:1:length(equip)
                k = equip(j);
                if strcmp(networkNames{net},'Electrical') 
                    if strcmp(Plant.Generator(k).Source,'Renewable')% subtract renewable generation 
                        QP.beq(req) = QP.beq(req) - Forecast.Renewable(t,k); %put renewable generation into energy balance at correct node
                    end
                end
                if ~isempty(FirstProfile) && ismember(Plant.Generator(k).Type,{'Electric Storage';'Thermal Storage';})
                    if isfield(Plant.Generator(k).QPform.output,out)
                        loss = (Plant.Generator(k).QPform.Stor.SelfDischarge*Plant.Generator(k).QPform.Stor.UsableSize*Plant.Generator(k).QPform.Stor.DischEff);
                        d_SOC = FirstProfile(t+1,k) - IC(k);%positive means charging
                        QP.beq(req) = QP.beq(req) + (d_SOC/dt(t)+loss)*Plant.Generator(k).QPform.Stor.DischEff;
                        QP.b(QP.Organize.Inequalities(k)) = -(1/Plant.Generator(k).QPform.Stor.ChargeEff-Plant.Generator(k).QPform.Stor.DischEff)*(d_SOC/dt(t)+loss);
                    end
                end
                if isfield(Plant.Generator(k).QPform,'constDemand') && isfield(Plant.Generator(k).QPform.constDemand,out)
                    QP.constDemand.(out).req(n,1) = req;
                    QP.constDemand.(out).load(n,k) = Plant.Generator(k).QPform.constDemand.(out);
                end
            end
        end
    end
end  

for i = 1:1:nB%Building inequalities & equality
    states = QP.Organize.States{nG+nL+i};
    req = QP.Organize.Building.Electrical.req(i);
    QP.beq(req) = QP.beq(req) + Forecast.Building.E0(t,i);%Equipment and nominal Fan Power
    %Heating Equality
    req = QP.Organize.Building.DistrictHeat.req(i);
    QP.Organize.Building.H_Offset(1,i) = Plant.Building(i).QPform.UA*(Forecast.Building.Tzone(t,i)-Forecast.Building.Tset_H(t,i)) - Forecast.Building.H0(t,i);
    QP.beq(req) = QP.beq(req) - QP.Organize.Building.H_Offset(1,i);%Heating load
    QP.lb(states(2)) = QP.lb(states(2)) + QP.Organize.Building.H_Offset(1,i);
    
    %Cooling Equality
    req = QP.Organize.Building.DistrictCool.req(i);
    QP.Organize.Building.C_Offset(1,i) = Plant.Building(i).QPform.UA*(Forecast.Building.Tset_C(t,i)-Forecast.Building.Tzone(t,i)) - Forecast.Building.C0(t,i);
    QP.beq(req) = QP.beq(req) - QP.Organize.Building.C_Offset(1,i);%Cooling Load
    QP.lb(states(3)) = QP.lb(states(3)) + QP.Organize.Building.C_Offset(1,i);
    
    %heating inequality
    QP.b(QP.Organize.Building.r(i)) = Plant.Building(i).QPform.UA*Forecast.Building.Tset_H(t,i) + Plant.Building(i).QPform.Cap*T(1,i)/(3600*dt(t));
    QP.A(QP.Organize.Building.r(i),states(1)) = (Plant.Building(i).QPform.UA+Plant.Building(i).QPform.Cap/(3600*dt(t)));
    %cooling inequality
    QP.b(QP.Organize.Building.r(i)+1) = -Plant.Building(i).QPform.UA*Forecast.Building.Tset_C(t,i) + Plant.Building(i).QPform.Cap*T(1,i)/(3600*dt(t));
    QP.A(QP.Organize.Building.r(i)+1,states(1)) = -Plant.Building(i).QPform.UA - Plant.Building(i).QPform.Cap/(3600*dt(t));
    %penalty states (excess and under temperature
    QP.b(QP.Organize.Building.r(i)+2) = Forecast.Building.Tmax(t,i);%upper buffer inequality, the longer the time step the larger the penaly cost on the state is.
    QP.b(QP.Organize.Building.r(i)+3) = -Forecast.Building.Tmin(:,i);%lower buffer inequality
end

%% Update upper and lower bounds based on ramping constraint (if applicable)
if ~isempty(FirstProfile)
    for i = 1:1:nG
        states = QP.Organize.States{i};
        if ismember(Plant.Generator(i).Type,{'Utility';'Electric Generator';'CHP Generator';'Chiller';'Heater';})%all generators and utilities
            if MinPower(i)>sum(QP.lb(states))
                QP.Organize.Dispatchable(i) = 0; %can't shut off
                %raise the lower bounds
                Padd = MinPower(i) - sum(QP.lb(states));
                for f = 1:1:length(states) %starting with lb of 1st state in generator, then moving on
                    gap = QP.ub(states(f)) - QP.lb(states(f));
                    if gap>0
                        gap = min(gap,Padd);
                        QP.lb(states(f)) = QP.lb(states(f)) + gap;
                        Padd = Padd-gap;
                    end
                end
            end
            if MaxPower(i)<Plant.Generator(i).Size %lower the upper bounds
                Psub = sum(QP.ub(states)) - MaxPower(i);
                for f = length(states):-1:1 %starting with lb of 1st state in generator, then moving on
                    if QP.ub(states(f))>0
                        gap = min(QP.ub(states(f)),Psub);
                        QP.ub(states(f)) = QP.ub(states(f)) - gap;
                        Psub = Psub-gap;
                    end
                end
            end
        elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
            %update storage output range to account for what is already scheduled
            if isfield(Plant.Generator(i).QPform.output,'H')
                req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            elseif isfield(Plant.Generator(i).QPform.output,'E')
                req = QP.Organize.Balance.Electrical;
            elseif isfield(Plant.Generator(i).QPform.output,'C')
                req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
            end
            s = states(1);
            QP.lb(s) = QP.lb(s) - StorPower(i);%change in storage for this power output
            QP.ub(s) = QP.ub(s) - StorPower(i);%change in storage for this power output
            chargingSpace = max(0,sum((Plant.Generator(i).QPform.Stor.UsableSize-FirstProfile(t+1,i)).*QP.Aeq(req,s))/dt(t));%you can't charge more than you have space for
            QP.lb(s) = max(QP.lb(s), -chargingSpace);
            QP.ub(s) = min(QP.ub(s), FirstProfile(t+1,i)/dt(t));%you cant discharge more than you have stored            
            if ismember(Plant.Generator(i).Type,{'Hydro Storage';})
                %update the inequality with (Power change from nominal)*conversion -Outflow above nominal < nominal spill flow
            end
            if QP.Organize.SpinRow(i)~=0 %update spinning reserve (max additional power from storage
                QP.beq(QP.Organize.SpinRow{i}) = QP.ub(s);
            end
        elseif ismember(Plant.Generator(i).Type,{'Hydro Storage';})
                %update the inequality with Power*conversion - Outflow above nominal < nominal OutFlow
                %change lb to 0 and ub to max gen
        end
        QP.lb(states) = min(QP.lb(states),QP.ub(states)); %fix rounding errors
    end
end

%% Adding cost for Line Transfer Penalties to differentiate from spillway flow
if ismember('Hydro',networkNames) && ~isempty(Plant.subNet.Electrical.lineNames)
    for k = 1:1:length(Plant.subNet.Electrical.lineNames)
        line = Plant.subNet.Electrical.lineNumber(k);
        states = QP.Organize.States{nG+line};
        if length(states)>1 %bi-directional transfer with penalty states
            s2 = states(2):QP.Organize.t1States:(states(2));%penalty from a to b
            s3 = states(3):QP.Organize.t1States:(states(3));% penalty from b to a
            QP.f(s2) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
            QP.f(s3) = 0.001; %Adding .1 cent/kWhr cost to every Electric line penalty to differentiate from spill flow
        end
    end 
end 

%% update costs
H = diag(QP.H);%convert to vector
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i}; %states of this generator
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') 
            QP.f(states(1)) = QP.f(states(1))*scaleCost(i)*dt(t);
            if Plant.Generator(i).VariableStruct.SellBackRate == -1
                QP.f(states(2)) = QP.f(states(2))*scaleCost(i)*dt(t); %sellback is a fixed percent of purchase costs (percent set wehn building matrices)
            else
                %constant sellback rate was taken care of when building the matrices
            end
        elseif ismember(Plant.Generator(i).Type,{'Electric Generator';'CHP Generator';'Utility';'Chiller';'Heater';})%all generators and utilities
            H(states) = H(states)*scaleCost(i)*dt(t);
            QP.f(states) = QP.f(states)*scaleCost(i)*dt(t);
        elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})%all  storage
            if strcmp(Plant.Generator(i).Source,'Electricity')
                type = 'E';
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                type = 'H';
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                type = 'C';
            elseif strcmp(Plant.Generator(i).Source,'Water')
                type = 'W';
            end
            %penalize deviations from the expected storage output
            %don't penalize spinning reserve state
            if isempty(FirstProfile)
                QP.f(states(1)) = marginal.(type);
                if any(strcmp(networkNames,'Electrical')) && strcmp(type,'H') % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                    H(states(1)) = 0;
                elseif strcmp(type,'W') %Hydro needs to have high cost, but not too high that it never uses it; set arbitrary number
                    H(states(1)) = 1e5*marginal.(type);
                else 
                    H(states(1)) = 1e8*marginal.(type);
                end
            else
                a = 4; %sets severity of quadratic penalty
                MaxCharge = Plant.Generator(i).QPform.Ramp.b(1); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                if Plant.Generator(i).QPform.Stor.UsableSize>FirstProfile(t+1,i)
                    MaxCharge = min((Plant.Generator(i).QPform.Stor.UsableSize-FirstProfile(t+1,i))/dt(t),MaxCharge);
                end
                QP.f(states(1)) = marginal.(type);
                if MaxCharge == 0
                    H(states(1)) = 0;
                else
                    H(states(1)) = a*2*QP.f(states(1))/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
                QP.f(states(2)) = 1e-6; %small penalty on charging state so that if there are no other generators it keeps the charging constraint as an equality
            end
        end
    end
end 
QP.constCost = QP.constCost.*scaleCost*dt(t);
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(1,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    SRancillary = QP.Organize.SpinReserveStates(1,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt(t)./(Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt(t)./(0.05*sum(Forecast.Demand.E(t,:),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt(t); % $0.05 per kWh
end
QP.H = diag(H);%convert back to matrix