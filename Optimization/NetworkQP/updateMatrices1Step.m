function QP = updateMatrices1Step(QP,Forecast,Renewable,scaleCost,dt,EC,StorPower,MinPower,MaxPower,T)
% update the equalities with the correct demand, and scale fuel and electric costs
% EC is the expected end condition at this time stamp (can be empty)
% StorPower is the expected output/input of any energy storage at this timestep (can be empty)
% MinPower and MaxPower define the range of this generator at this timestep
%T is the building temperatures
global Plant
QP.solver = Plant.optimoptions.solver;
nG = length(Plant.Generator);
nB = length(Plant.Building);
nL = length(QP.Organize.States)-nB-nG;
marginal = instantMarginalCost(EC,scaleCost);%update marginal cost
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
                QP.b(QP.Organize.SpinReserve) = -Plant.optimoptions.SpinReservePerc/100*sum(Forecast.E,2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
            end
        elseif strcmp(networkNames{net},'DistrictHeat')
            out = 'H';
        elseif strcmp(networkNames{net},'DistrictCool')
            out = 'C';
        end
        for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
            equip = Plant.subNet.(networkNames{net}).Equipment{n};
            req = QP.Organize.Balance.(networkNames{net})(n);
            loads = nonzeros(QP.Organize.Demand.(networkNames{net}){n});
            if ~isempty(loads)%need this to avoid the next line when there is no load
                QP.beq(req) = sum(Forecast.(out)(:,loads),2); %multiple demands can be at the same node
            end
            for j = 1:1:length(equip)
                i = equip(j);
                if strcmp(networkNames{net},'Electrical') && strcmp(Plant.Generator(i).Source,'Renewable')% subtract renewable generation 
                    QP.beq(req) = QP.beq(req) - Renewable(i); %put renewable generation into energy balance at correct node
                end
                if ~isempty(StorPower) && ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
                    if isfield(Plant.Generator(i).QPform.output,out)
                        QP.beq(req) = QP.beq(req) - StorPower(i); %%remove expected storage output from beq
                    end
                end
            end
        end
    end
end  

for i = 1:1:nB%Building inequalities & equality
    states = QP.Organize.States{nG+nL+i};
    req = QP.Organize.Building.req(i);
    QP.beq(req) = QP.beq(req) + Forecast.Building(i).E0 - Plant.Building(i).QPform.H2E*Forecast.Building(i).H0 - Plant.Building(i).QPform.C2E*Forecast.Building(i).C0;%electric equality
    
    QP.b(QP.Organize.Building.r(i)) = Plant.Building(i).QPform.UA*Forecast.Building(i).Tset_H - Forecast.Building(i).H0 + Plant.Building(i).QPform.Cap*T(i)/(3600*dt);%heating inequality
    QP.A(QP.Organize.Building.r(i),states(1)) = (Plant.Building(i).QPform.UA+Plant.Building(i).QPform.Cap/(3600*dt));% Heating>= H0 + UA*(Ti-Tset) + Cap*(Ti - T(i-1))/dt where dt is in seconds
    QP.b(QP.Organize.Building.r(i)+1) = -Plant.Building(i).QPform.UA*Forecast.Building(i).Tset_C - Forecast.Building(i).C0 - + Plant.Building(i).QPform.Cap*T(i)/(3600*dt);%cooling inequality
    QP.A(QP.Organize.Building.r(i)+1,states(1)) = -(Plant.Building(i).QPform.UA+Plant.Building(i).QPform.Cap/(3600*dt));% Cooling>= C0 + UA*(Tset-Ti) + Cap*(T(i-1)-Ti)/dt where dt is in seconds
    
    QP.b(QP.Organize.Building.r(i)+2) = (Forecast.Building(i).Tset + Plant.Building(i).VariableStruct.Comfort/2);%upper buffer inequality, the longer the time step the larger the penaly cost on the state is.
    QP.b(QP.Organize.Building.r(i)+3) = -(Forecast.Building(i).Tset - Plant.Building(i).VariableStruct.Comfort/2);%lower buffer inequality
end

%% Update upper and lower bounds based on ramping constraint (if applicable)
if ~isempty(MinPower)
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
        elseif ~isempty(StorPower) && ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydro Storage';})
            %update storage output range to account for what is already scheduled
            QP.lb(states) = QP.lb(states) - StorPower(i);
            QP.ub(states) = QP.ub(states) - StorPower(i);            
            if ismember(Plant.Generator(i).Type,{'Hydro Storage';})
                %update the inequality with (Power change from nominal)*conversion -Outflow above nominal < nominal spill flow
            end
            if QP.Organize.SpinRow(i)~=0 %update spinning reserve (max additional power from storage
                QP.beq(QP.Organize.SpinRow{i}) = sum(QP.ub(states(1:end-1))) - StorPower(i);
            end
        elseif ismember(Plant.Generator(i).Type,{'Hydro Storage';})
                %update the inequality with Power*conversion - Outflow above nominal < nominal OutFlow
                %change lb to 0 and ub to max gen
        end
    end
end

%% update costs
H = diag(QP.H);%convert to vector
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i}; %states of this generator
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') && Plant.Generator(i).VariableStruct.SellBackRate>0
            QP.f(states(1)) = QP.f(states(1))*scaleCost(i)*dt;
            QP.f(states(2)) = min(0.9999*QP.f(states(1)),Plant.Generator(i).VariableStruct.SellBackRate*dt); %make sure sellback rate is less than purchase rate
        elseif ismember(Plant.Generator(i).Type,{'Electric Generator';'CHP Generator';'Utility';'Chiller';'Heater';})%all generators and utilities
            H(states) = H(states)*scaleCost(i)*dt;
            QP.f(states) = QP.f(states)*scaleCost(i)*dt;
        else
            if strcmp(Plant.Generator(i).Source,'Electricity')
                type = 'Electrical';
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                type = 'DistrictHeat';
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                type = 'DistrictCool';
            elseif strcmp(Plant.Generator(i).Source,'Water')
                type = 'Hydro';
            end
            %penalize deviations from the expected storage output
            %don't penalize spinning reserve state
            if isempty(EC)
                QP.f(states(1)) = marginal.(type);
                if any(strcmp(networkNames,'Electrical')) && strcmp(type,'DistrictHeat') % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                    H(states(1)) = 0;
                else
                    H(states(1)) = 1e8*marginal.(type);
                end
            else
                a = 4; %sets severity of quadratic penalty
                MaxCharge =min((Plant.Generator(i).QPform.Stor.UsableSize-EC(i))/dt,Plant.Generator(i).QPform.Ramp.b(1)); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                QP.f(states(1)) = marginal.(type);
                if MaxCharge == 0
                    H(states(1)) = 0;
                else
                    H(states(1)) = a*2*QP.f(states(1))/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
            end
        end
    end
end 
QP.constCost = QP.constCost.*scaleCost*dt;
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(1,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    SRancillary = QP.Organize.SpinReserveStates(1,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt./(Plant.optimoptions.SpinReservePerc/100*sum(Forecast.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt./(0.05*sum(Forecast.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);%convert back to matrix