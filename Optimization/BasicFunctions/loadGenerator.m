function loadGenerator% Loads generators for economic dispatch
%% this function identifies the values that will be used to represent each generator in the quadratic optimizations
global Plant
nG = length(Plant.Generator);
for i = 1:1:nG
    Plant.Generator(i).QPform = {}; %delete this line when you begin using the gui again
    if isempty(Plant.Generator(i).QPform)%only load generators that have not been loaded yet. New run, new generator, or edited generator
        typeNoSpace = char(Plant.Generator(i).Type(~isspace(char(Plant.Generator(i).Type))));
        [Plant.Generator(i).QPform, Plant.Generator(i).VariableStruct.dX_dt,SS] = eval(strcat('load',typeNoSpace,'(Plant.Generator(i))'));
        if ~isempty(SS)
            SSi(i) = SS;
        end
    end
end 
if ~exist('SSi','var')
    SSi = [];
end
agregateSSmodel(SSi)
end%Ends function loadGenerator

function [QPform, dX_dt,SSi] = loadUtility(Gen)
util = Gen.VariableStruct;
dX_dt = inf;
SSi =[];
if any(strcmp(Gen.Source, {'NG'}))%if it is a fuel utility
    QPform = [];
    QPform.output = [];
    QPform.states = [];
else
    QPform.states = {'X'};
    QPform.X.H = 0;
    QPform.X.f = 1;
    QPform.X.lb = util.MinImportThresh;
    QPform.X.ub = inf;% no sellback allowed (only 1 state)
   
    if util.MinImportThresh<=0 && (util.SellBackPerc>0 || util.SellBackRate>0) %add sell back state
        QPform.states = {'X';'Y'};
        if util.SellBackPerc>0
            QPform.Y.f = -min(util.SellBackPerc/100,1-1e-6);%ensure less than 1, so no issues with pass through power
        elseif util.SellBackRate>0
            QPform.Y.f = -1;%constant sell back rate
        end
        QPform.Y.H = 0;
        QPform.Y.lb = 0;
        QPform.Y.ub = inf;
        if strcmp(Gen.Source, 'Electricity')
            QPform.output.E = [1;-1];
        elseif strcmp(Gen.Source, 'Heat')%loads the parameters for a distric heating supply. 
            QPform.output.H = [1;-1];
        elseif strcmp(Gen.Source, 'Cooling')%loads the parameters for a distric cooling supply. 
            QPform.output.C = [1;-1];
        else 
        end
    else
        if strcmp(Gen.Source, 'Electricity')
            QPform.output.E = 1;
        elseif strcmp(Gen.Source, 'Heat')%loads the parameters for a distric heating supply. 
            QPform.output.H = 1;
        elseif strcmp(Gen.Source, 'Cooling')%loads the parameters for a distric cooling supply. 
            QPform.output.C = 1;
        else 
        end
    end
    
end
end%Ends function loadUtility

function [QPform, dX_dt,SSi] = loadElectricGenerator(Gen)
% this function loads the parameters for an electric generator generators
[QPform,dX_dt,SSi] = loadCHPGenerator(Gen);
end%Ends function loadElectricGenerator

function [QPform,dX_dt,SSi] = loadChiller(Gen)
global Plant
if Plant.optimoptions.sequential
    % this function loads the parameters for a chiller where the cost is in
    % kW of electric power. It must be optimized first and the resulting
    % demand added to the electric load
    [QPform,dX_dt,SSi] = loadCHPGenerator(Gen);
    %% can't do this method with an absorption chiller in the mix
else
    %this loads a segmented convex fit of an electric chiller, where the eletric load 
    %shows up in the electric energy balance
        [QPform,dX_dt,SSi] = segmentedChiller(Gen);
%         [QPform,dX_dt,SSi] = constantChiller(Gen);
end
end%Ends function loadChiller

function [QPform,dX_dt,SSi] = constantChiller(Gen)
%loads a convex segmented COP fit of an electric or absorption chiller
global Plant
if strcmp(Gen.Source,'Electricity')
    source = 'E';
else
    source = 'H';
end
maxCOP = max(Gen.Output.Cooling);
LB = Gen.VariableStruct.Startup.Cooling(end); %chiller
UB = Gen.Size;
[dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,[]);
dX_dt = dX_dt/Plant.optimoptions.scaletime;  
QPform.Ramp.b = [dX_dt;dX_dt]; %-output1+output2=ramp up %output1-output2=-rampdown
QPform.states(1,1) = {'A'};
QPform.states(1,2) = {'A'};
QPform.A.H = 0;
QPform.A.f = 0;
QPform.A.lb = 0;
QPform.A.ub = UB;
QPform.output.C = 1;
QPform.output.(source) = -1/maxCOP;

%% Fit B has a constant electric demand
xi = min(length(Gen.Output.Capacity)-2,nnz(Gen.Output.Capacity<0.5*Gen.Output.Capacity(end))); %fit the top 50% of output
FitB = polyfit(Gen.Output.Capacity(xi:end)*UB,Gen.Output.Capacity(xi:end)*UB./Gen.Output.Cooling(xi:end),1);%calculating fit of input vs output
QPform.A.H(2) = 0;
QPform.A.f(2) = 0;
QPform.A.lb(2) = LB;
QPform.A.ub(2) = UB;
QPform.output.C(2) = 1;
QPform.constDemand.(source) = FitB(2);%constant electrical demand (kWe) when using fitB
QPform.output.(source)(2) = -FitB(1);%marginal kWe for each kWthermal
end%Ends function constantChiller

function [QPform,dX_dt,SSi] = segmentedChiller(Gen)
%loads a convex segmented COP fit of an electric or absorption chiller
global Plant
if strcmp(Gen.Source,'Electricity')
    source = 'E';
else
    source = 'H';
end
[~,index] = max(Gen.Output.Cooling);
LB = Gen.VariableStruct.Startup.Cooling(end); %chiller
UB = Gen.Size;
n = min(3,length(Gen.Output.Cooling)-index+1); %number of linear segments for the chiller COP fit
[dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,[]);
dX_dt = dX_dt/Plant.optimoptions.scaletime;  
QPform.Ramp.b = [dX_dt;dX_dt]; %-output1+output2=ramp up %output1-output2=-rampdown
letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';};
QPform.states(1:n,1) = letters(1:n);
QPform.states(1:n,2) = letters(1:n);
xi = 2;
ub = max(Gen.Output.Capacity(index)*UB,LB + (UB-LB)/n);
for i = 1:1:n
    QPform.(letters{i}).H = 0;
    QPform.(letters{i}).f = 0;
    QPform.(letters{i}).lb = 0;
    if i ==1
        QPform.(letters{i}).ub = ub;
    else
        ub = ub + (UB-QPform.(letters{1}).ub(1))/(n-1);
        QPform.(letters{i}).ub = (UB-QPform.(letters{1}).ub(1))/(n-1);
    end

    QPform.(letters{i}).H(2) = 0;
    QPform.(letters{i}).f(2) = 0;
    if i ==1
        QPform.(letters{i}).lb(2) = LB;
        QPform.(letters{i}).ub(2) = QPform.(letters{i}).ub;
    else
        QPform.(letters{i}).lb(2) = 0;
        QPform.(letters{i}).ub(2) = (UB-QPform.(letters{1}).ub(1))/(n-1);
    end
    xend = xi;
    while Gen.Output.Capacity(xend)*UB<ub && xend<length(Gen.Output.Capacity)
        xend = xend+1;
    end
    COPavg = mean(interp1(Gen.Output.Capacity(xi:xend)*UB,Gen.Output.Cooling(xi:xend),linspace(max(Gen.Output.Capacity(xi)*UB,ub-QPform.(letters{i}).ub(1)),ub,10)'));
    FitB = polyfit(Gen.Output.Capacity(xi:xend)*UB,Gen.Output.Capacity(xi:xend)*UB./Gen.Output.Cooling(xi:xend),1);%calculating fit of input vs output
    if i ==1
        QPform.constDemand.(source) = FitB(2);%constant  demand (kW) 
    end
    QPform.output.C(i,1:2) = 1;
    QPform.output.(source)(i,1) = -1/COPavg;
    QPform.output.(source)(i,2) = -FitB(1);
    xi = max(1,xend-1);
    
%     if i ==1
%         figure(1)
%         plot(Gen.Output.Capacity*UB,Gen.Output.Capacity*UB./Gen.Output.Cooling)
%     elseif i ==2
%         hold on
%         plot(Gen.Output.Capacity(xi:xend)*UB, Gen.Output.Capacity(xi:xend)*UB*FitB(1)+FitB(2),'r')
%     elseif i == 3
%         plot(Gen.Output.Capacity(xi:xend)*UB, Gen.Output.Capacity(xi:xend)*UB*FitB(1)+FitB(2),'g')
%     elseif i ==4
%         plot(Gen.Output.Capacity(xi:xend)*UB, Gen.Output.Capacity(xi:xend)*UB*FitB(1)+FitB(2),'k')
%     end
end
end%Ends function segmentedChiller

function [QPform,dX_dt,SSi] = loadCHPGenerator(Gen)
% this function loads the parameters for a combined heat and power
% generator, regular electric generator, or chiller
global Plant
UB = Gen.Size;
if isfield(Gen.Output,'Cooling')&& Gen.Output.Cooling(end)>0
    LB = Gen.VariableStruct.Startup.Cooling(1); %chiller
    costTerms = GenCosts(Gen,LB,UB,'Cooling');
    QPform.output.C = 1;
else
    LB = Gen.VariableStruct.Startup.Electricity(end); %electric or CHP generator 
    costTerms = GenCosts(Gen,LB,UB,'Electricity');
    QPform.output.E = 1;
end
if isfield(Gen.Output,'Heat')&& Gen.Output.Heat(end)>0
    Hratio = Gen.Output.Heat(end)/Gen.Output.Electricity(end);
    QPform.output.H = Hratio;
else Hratio = [];
end
if license('test','Control_Toolbox')
    [dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,Hratio);
else dX_dt = Gen.Size/4; SSi = []; %assume 4 hours from off to peak
end
dX_dt = dX_dt/Plant.optimoptions.scaletime;  
QPform.Ramp.b = [dX_dt;dX_dt]; %-output1+output2=ramp up %output1-output2=-rampdown
QPform.constCost = costTerms.Intercept(4);
if costTerms.P == UB %linear fit use 1 state
    QPform.states = {'X'};
    QPform.X.H = 0;
    QPform.X.f = costTerms.Convex(1);
    QPform.X.lb = 0;%relax lower bound until have figured out what will be off or on
    QPform.X.ub = UB;

else %piecewise quadratic fit
    QPform.states = {'Y';'Z';};%beta, gamma
    
    QPform.Y.H = 0;
    QPform.Y.f = costTerms.Convex(1);
    QPform.Y.lb = 0;%relax lower bound until have figured out what will be off or on
    QPform.Y.ub = costTerms.P;

    QPform.Z.H = 2*costTerms.Convex(3);
    QPform.Z.f = costTerms.Convex(2);
    QPform.Z.lb = 0;
    QPform.Z.ub = UB-costTerms.P;
end
%% FitB
if costTerms.I == UB %linear fit use 1 state
    QPform.states(1,2) = {'X'};
    QPform.X.H(2) = 0;
    QPform.X.f(2) = costTerms.Intercept(1);
    QPform.X.lb(2) = LB;%relax lower bound until have figured out what will be off or on
    QPform.X.ub(2) = UB;
else
    QPform.states(1:2,2) = {'Y';'Z'};%x, beta, gamma
    QPform.Y.H(2) = 0;
    QPform.Y.f(2) = costTerms.Intercept(1);
    QPform.Y.lb(2) = LB;
    QPform.Y.ub(2) = costTerms.I;

    QPform.Z.H(2) = 2*costTerms.Intercept(3);
    QPform.Z.f(2) = costTerms.Intercept(2);
    QPform.Z.lb(2) = 0;
    QPform.Z.ub(2) = UB-costTerms.I;  
end
end%Ends function loadCHPGenerator

function [QPform,dX_dt,SSi] = loadHeater(Gen)
% this function loads the parameters for a heater generator.
global Plant
LB = Gen.VariableStruct.Startup.Heat(end);
[dX_dt, SSi] = RampRateCalc(Gen.VariableStruct.StateSpace, LB, Gen.Size, {});
dX_dt = dX_dt/Plant.optimoptions.scaletime;
capacity = Gen.Output.Capacity*Gen.Size;
efficiency = Gen.Output.Heat;
costLinear = 1/mean(efficiency(capacity>=LB));%efficiency term  put in cost (will later be scaled by utility cost.

QPform.states = {'X'};
QPform.output.H = 1;
QPform.Ramp.b = [dX_dt; dX_dt];%ramp up, ramp down

QPform.X.H = 0;
QPform.X.f = costLinear;
QPform.X.lb = 0;
QPform.X.ub = Gen.Size;

QPform.states(1,2) = {'X'};
QPform.X.H(2) = 0;%fit B
QPform.X.f(2) = costLinear;%fit B
QPform.X.lb(2) = LB; %fit B
QPform.X.ub(2) = Gen.Size;%fit B
end%Ends function loadHeater

function [QPform,dX_dt,SSi] = loadElectricStorage(Gen)
%this function loads all the parameters needed for electric storage
[QPform,dX_dt,SSi] = Storage(Gen);
end%Ends function loadElectricStorage

function [QPform,dX_dt,SSi] = loadThermalStorage(Gen)
%loads all types of thermal storage
[QPform,dX_dt,SSi] = Storage(Gen);
end%Ends function loadThermalStorage

function [QPform,dX_dt,SSi] = Storage(Gen)
%this function just directs to either hot or cold thermal storage
%if we can get rid of the CS, HS structures then this function can be
%eliminated, because all storage can be handled the same way.
global Plant
SSi =[];
Stor.Size = Gen.Size*Plant.optimoptions.scaletime;
Stor.SelfDischarge  = Gen.VariableStruct.SelfDischarge;% SelfDischarge per hour (fraction of total charge)
if isfield(Gen.VariableStruct, 'EnStoreType') %if its thermal storage
    Stor.PeakDisch = (Gen.VariableStruct.DischRatePerc/100*Gen.Size); %Thermal kW out
    Stor.PeakCharge = (Gen.VariableStruct.FillRatePerc/100*Gen.Size); %Thermal kW in
    Stor.ChargeEff = Gen.VariableStruct.ChargeEff;
    Stor.DischEff = Gen.VariableStruct.DischargeEff;
    Stor.UsableSize  = Stor.Size; % usable size 
    if strcmp(Gen.VariableStruct.EnStoreType, 'ColdTES')
        QPform.output.C = 1; 
    elseif strcmp(Gen.VariableStruct.EnStoreType, 'HotTES')
        QPform.output.H = 1;
        Stor.ChargeEff = 1; %set to ideal 
        Stor.DischEff = 1;
    end
    
else %electric battery
    QPform.output.E = 1; 
    Stor.Voltage = Gen.VariableStruct.Voltage;
    DischCurrent = Gen.VariableStruct.PeakDisch.*Gen.Size/Stor.Voltage*1000;
    Stor.DischResistScaled = (100/DischCurrent)*Gen.VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    DischVoltLoss = DischCurrent*Stor.DischResistScaled; %keep in mind when calculating loss as function of discharge current
    ChargeCurrent = Gen.VariableStruct.PeakCharge*Gen.Size/Stor.Voltage*1000;
    Stor.ChargeResistScaled = (100/ChargeCurrent)*Gen.VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    ChargeVoltLoss = ChargeCurrent*Stor.ChargeResistScaled;
    Stor.ChargeEff = Stor.Voltage/(Stor.Voltage+ChargeVoltLoss);
    Stor.DischEff = (Stor.Voltage-DischVoltLoss)/Stor.Voltage;
    Stor.PeakDisch = (DischCurrent*Stor.Voltage*Stor.DischEff/1000); %PeakDischPower kW out
    Stor.PeakCharge = (ChargeCurrent*Stor.Voltage/Stor.ChargeEff/1000); % PeakChargePower kW in
    Stor.UsableSize  = Stor.Size*(Gen.VariableStruct.MaxDOD/100); % usable size (state x is 0_> usable size, must calculate this initial condition each time from voltage & SOC
end

QPform.Stor = Stor;
dX_dt = Stor.PeakDisch;% storage discharge constraint in kW

QPform.states = {'X';}; %state of charge, charging power, no buffers
QPform.X.lb = 0;
QPform.X.ub = Stor.UsableSize;
QPform.X.H = 0;
QPform.X.f = 0;
QPform.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
        
        
a = (1/Stor.ChargeEff - Stor.DischEff);
if a~=0 %not ideal storage, add charging state
   QPform.states = {'X'; 'Y'}; %state of charge, charging power, no buffers
    QPform.link.ineq = [a -1]; % (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
    QPform.link.bineq = 0;
    QPform.Y.lb = 0;
    QPform.Y.ub = Stor.PeakCharge;%the limit on how much charging power can be delivered is handled by the generators' limits, so put inf here to prevent redundancy
    QPform.Y.H = 0;%the cost of the charging power is handled by the generators 
    QPform.Y.f = 0;
end    
if Plant.optimoptions.Buffer ~= 0 %buffer states
    if a==0 %ideal storage, ignore charging state
        QPform.states = {'X';'U';'L'};%SOC(t+1), charging power, upper buffer, lower buffer
        QPform.link.ineq = [-1 0 -1; 1 -1 0];%-SOC(t)-lowerbuffer<-0.2  and %SOC-upperbuffer<0.8
        QPform.link.bineq = [0;0;]; %% note: the magnitude of the buffer is set later in findBuffer
    else
        QPform.states = {'X';'Y';'U';'L'};%SOC(t+1), charging power, upper buffer, lower buffer
        QPform.link.ineq = [a, -1, 0, 0];% (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
        QPform.link.ineq = [QPform.link.ineq; -1 0 0 -1];%-SOC(t)-lowerbuffer<-0.2
        QPform.link.ineq = [QPform.link.ineq; 1 0 -1 0];%SOC-upperbuffer<0.8
        QPform.link.bineq = [0;0;0;];
    end
    QPform.U.lb = 0;
    QPform.U.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    QPform.U.H = 0;
    QPform.U.f = 0;

    QPform.L.lb = 0;
    QPform.L.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    QPform.L.H = 0;
    QPform.L.f = 0;
end
end%Ends function Storage

function [QPform,dX_dt,SSi] = loadSolar(Gen)
%PV solar
QPform.output.E = 1;%there are no states or outputs for solar because renewable outputs are handled on the demand side
QPform.states = [];
dX_dt = inf;
SSi = [];
end%Ends function loadSolar

function [QPform, dX_dt,SSi] = loadHydroStorage(Gen)
% this function loads the parameters for a hydroelectric plant.
global Plant
SSi =[];
Eff = Gen.VariableStruct.MaxGenCapacity/(Gen.VariableStruct.MaxGenFlow*Gen.VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW
QPform.Stor.Size = Gen.Size*Plant.optimoptions.scaletime;
QPform.Stor.SelfDischarge  = 0; %needs to be evaporative losses
%Put this at 1/2 full capacity until the MaxHead-MinHead/MaxHead is fixed
%Currently, with the previous solution, some plants lost a vast majority of usable space
QPform.Stor.UsableSize  = QPform.Stor.Size*0.5;%*((Gen.VariableStruct.MaxHead-Gen.VariableStruct.MinHead)/Gen.VariableStruct.MaxHead);
QPform.Stor.Power2Flow = 1/(Eff*Gen.VariableStruct.MaxHead*84.674);%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 84.674 kJ/ (1000ft^3*ft)
QPform.output.E = [1 0];
QPform.output.W = 0;

QPform.states = {'X';'Y'}; %Power and state of charge
QPform.X.lb = 0;
QPform.X.ub = Gen.VariableStruct.MaxGenCapacity;
QPform.X.H = 0;
QPform.X.f = 0;

QPform.Y.lb = 0;
QPform.Y.ub = QPform.Stor.UsableSize;
QPform.Y.H = 0;
QPform.Y.f = 0;

QPform.Ramp.b = [Gen.VariableStruct.RampDown; Gen.VariableStruct.RampUp;];%change in power generation

QPform.link.eq = [QPform.Stor.Power2Flow, 0];%convert power to flow rate (1000 cfs)
QPform.link.beq = 0;
    
if Gen.VariableStruct.MaxSpillFlow>0
    QPform.states = {'X';'Y';'S'}; %add spill flow state
    QPform.S.lb = 0;
    QPform.S.ub = Gen.VariableStruct.MaxSpillFlow;
    QPform.S.H = 0;
    QPform.S.f = 0;
    QPform.link.eq = [QPform.Stor.Power2Flow, 0, 1];
    QPform.link.beq = 0;
    QPform.output.E = [1 0 0];
end
if Plant.optimoptions.Buffer ~= 0 %buffer states
    if Gen.VariableStruct.MaxSpillFlow>0
        QPform.states = {'X';'Y';'S';'U';'L'};
        QPform.link.eq = [QPform.Stor.Power2Flow, 0, 1, 0, 0];
        QPform.link.ineq = [0 -1 0 0  -1; 0 1 0 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        QPform.output.E = [1 0 0 0 0];
    else
        QPform.states = {'X';'Y';'U';'L'};
        QPform.link.eq = [QPform.Stor.Power2Flow, 0, 0, 0];
        QPform.link.ineq = [0 -1 0  -1; 0 1 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        QPform.output.E = [1 0 0 0];
    end
    QPform.link.bineq = [0;0;];

    QPform.U.lb = 0;
    QPform.U.ub = 0;
    QPform.U.H = 0;
    QPform.U.f = 0;

    QPform.L.lb = 0;
    QPform.L.ub = 0;
    QPform.L.H = 0;
    QPform.L.f = 0;
end
dX_dt = Gen.VariableStruct.RampUp;
end%Ends function loadHydroStorage