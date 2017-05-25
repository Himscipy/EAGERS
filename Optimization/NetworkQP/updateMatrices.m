function QP = updateMatrices(QP,IC,Date,scaleCost,marginCost,Forecast,EC)
% QP is the set of optimization matrices 
% IC is the intial condition
% Time is the vector of time steps
% scale cost is the matrix multiplier of cost for each generator at each time step
% MarginCost is the marginal cost of generation (used in the final value of the energy storage)
% Demand is the forecasted loads
% Renewable is the forecasted uncontrollable generation
% EC is the end condition for the threshold optimization
global Plant
nG = length(Plant.Generator);
nS = length(Date)-1;
dt = (Date(2:end) - Date(1:end-1))*24;
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
nL = sum(nLinet);
%% update demands and storage self-discharge
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
        if Plant.optimoptions.SpinReserve
            QP.b(QP.Organize.SpinReserve) = -Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E,2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        end
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
    elseif strcmp(networkNames{net},'Hydro')
        out = 'W';
    end
    for i = 1:1:length(Plant.subNet.(networkNames{net})) %run through all the nodes in this network
        equip = Plant.subNet.(networkNames{net})(i).Equipment; %equipment at this node
        eq = QP.Organize.Balance.(networkNames{net})(i,:);%balance at this node
        load = QP.Organize.Demand.(networkNames{net})(i); %load at this node
        if load~=0
            QP.beq(eq) = sum(Forecast.Demand.(out)(:,load),2); %multiple demands can be at the same node, or none
        end
        for j = 1:1:length(equip)
            k = equip(j);
            if strcmp(networkNames{net},'Electrical') && any(QP.Renewable(:,k))% subtract renewable generation 
                QP.beq(eq) = QP.beq(eq) - QP.Renewable(:,k); %put renewable generation into energy balance at correct node
            end
            if ~isempty(strfind(Plant.Generator(k).Type,'Storage'))
                if isfield(Plant.Generator(k).OpMatA.output,out)
                    loss = dt*(Plant.Generator(k).OpMatA.Stor.SelfDischarge*Plant.Generator(k).OpMatA.Stor.UsableSize);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                end
            end
        end
    end
    if strcmp(networkNames{net},'Hydro')
        for i = 1:1:length(Plant.subNet.Hydro) %run through all the nodes in this network
            %% Node inflows (source/sink terms and upstream flow at time t-T ago if t<T)
            %need to be able to forecast sink/source
            I = find(strcmp(Plant.subNet.Hydro(i).nodes{1},Plant.Data.Hydro.Nodes));%column index of this node in the stored matrices of Data.Hydro.SourceSink and Data.Hydro.Inflow
            eq = QP.Organize.Balance.Hydro(i,:);%mass balance at this node
            QP.beq(eq)= QP.beq(eq) - Forecast.Hydro.SourceSink(I)./(12.1*dt); %source/sink term converted from 1000 ft^3/s to acre-ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
            K = Plant.subNet.Hydro(i).UpRiverSegments;%lines entering this node, i.e. this node is the downriver node
            for j = 1:1:length(K)
                Date2 = Date(2:end) - Plant.subNet.lineTime.Hydro(K(j))/24;
                n = nnz(Date2<Date(1));
                if n > 0 %if river segment is shorter than 1st time step, it doesn't need these constants, it uses initial condition & 1st step.
                    InFlow = getHydroFlows(Date2(1:n),K(j));%river segments flowing into this node from time before forecast.
                    QP.beq(eq(1:n))= QP.beq(eq(1:n)) - InFlow./(12.1*dt(1:n)); % Qupriver, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                end
            end
        end
    end
    nLcum = nLcum + nLinet(net);
end         

%update initial conditions
for i = 1:1:nG+nL
    if QP.Organize.IC(i)>0 
        nIC = nnz(QP.Organize.IC(1:i)); %order in IC
        QP.beq(nIC) = IC(i);
    end
end

%update costs
H = diag(QP.H);
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i};
        nt = length(states)/nS; %number of states per timestep
        if strcmp(Plant.Generator(i).Type,'Utility') && strcmp(Plant.Generator(i).Source,'Electricity') && Plant.Generator(i).VariableStruct.SellBackRate>0
            QP.f(states(1:2:end)) = QP.f(states(1:2:end)).*scaleCost(:,i).*dt;
            QP.f(states(2:2:end)) = QP.f(states(2:2:end)).*min(0.9999*scaleCost(:,i),Plant.Generator(i).VariableStruct.SellBackRate).*dt; %make sure sellback rate is less than purchase rate
        elseif isempty(strfind(Plant.Generator(i).Type,'Storage'))%all generators and utilities
            for t = 1:1:nS
                Xn = states((t-1)*nt+1:t*nt);
                H(Xn) = H(Xn)*scaleCost(t,i)*dt(t); 
                QP.f(Xn) = QP.f(Xn)*scaleCost(t,i)*dt(t); 
            end          
        else % storage systems
            if strcmp(Plant.Generator(i).Source,'Electricity')
                type = 'Electrical';
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                type = 'DistrictHeat';
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                type = 'DistrictCool';
            elseif strcmp(Plant.Generator(i).Source,'Water')
                type = 'Hydro';
            end
             %% update storage costs
            s_end = states(end)-nt+1; %final state of charge
            StorSize = Plant.Generator(i).OpMatA.X.ub;
            if isfield(Plant.Generator(i).OpMatA,'W') %has buffer
                BuffSize = Plant.Generator(i).OpMatA.W.ub;
            else BuffSize = 0;
            end
            if isempty(EC)
                if strcmp(Plant.Generator(i).Source,'Heat')
                    Max = 0.8*marginCost.DistrictHeat.Max;
                    Min = 0;
                elseif strcmp(Plant.Generator(i).Source,'Electricity')
                    Max = 1.25*marginCost.Electrical.Max;
                    Min = 0.95*marginCost.Electrical.Min;
                elseif strcmp(Plant.Generator(i).Source,'Cooling')
                    Max = 1.25*marginCost.DistrictCool.Max;
                    Min = 0.65*marginCost.DistrictCool.Min;
                elseif strcmp(Plant.Generator(i).Source,'Water')
                    Max = 1.25*marginCost.Hydro.Max;
                    Min = 0.65*marginCost.Hydro.Min;
                end
                a1 = -Max; % fitting C = a1*SOC + 0.5*a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                a2 = (Max - Min)/(StorSize);
                H(s_end) = a2;%quadratic final value term loaded into SOC(t=nS)  %quadprog its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = a1;%linear final value term loaded into SOC(t=nS)
                if BuffSize>0
                    for t = 1:1:nS
                        Xn = states(t*nt-1:t*nt);
                        QP.f(Xn) = Min;%this is the linear buffer term loaded into the lower & upper buffer
                        H(Xn) = (2*Max-Min)/BuffSize;%this is the quadratic buffer term loaded into the lower & upper buffer
                    end
                end
            else %used in Online loop when there is a target EC determined by dispatch loop
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
                rows = QP.Organize.Inequalities{i};
                nr = length(states)/nS;
                PeakChargePower = Plant.Generator(i).OpMatA.Ramp.b(1);
                dSOC_10perc = .1*PeakChargePower*(Date(end)-Date(1)); %energy (kWh) if charging at 10%
                H(s_end) = -2*marginCost.(type).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = -marginCost.(type).Min;%linear final value term loaded into SOC(t=nS)
                nIC = nnz(QP.Organize.IC(1:i)); %order in IC
                QP.beq(nIC) = IC(I)-EC(I);
                for t = 1:1:nS
                    Xn = states((t-1)*nt+1);
                    Rn = rows(t*nr);
                    QP.lb(Xn) = -EC(i);%change lb so that SOC = 0 coresponds to EC
                    QP.ub(Xn) = StorSize - EC(i);%change ub so that SOC = 0 coresponds to EC
                    QP.b(Rn-1) = -BuffSize + EC(i);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
                    QP.b(Rn) = StorSize-BuffSize - EC(i);%change lb so that SOC = 0 coresponds to EC
                end       
            end
        end
    end
end 
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(:,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    SRancillary = QP.Organize.SpinReserveStates(:,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt./(Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E,2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt./(0.05*sum(Forecast.Demand.E,2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);