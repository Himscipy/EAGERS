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
    if strcmp(Gen.Source, 'Electricity')
        QPform.output.E = 1;
    elseif strcmp(Gen.Source, 'Heat')%loads the parameters for a distric heating supply. 
        QPform.output.H = 1;
    elseif strcmp(Gen.Source, 'Cooling')%loads the parameters for a distric cooling supply. 
        QPform.output.C = 1;
    else 
    end
    if util.MinImportThresh<=0 && (util.SellBackPerc>0 || util.SellBackRate>0) %add sell back state
        QPform.states = {'X';'Y'};
        if util.SellBackPerc == 100 || util.SellBackRate>0
            QPform.Y.f = -1;%constant sell back rate
        else
            QPform.Y.f = -max(util.SellBackPerc/100,1-1e-6);%ensure less than 1, so no issues with pass through power
        end
        QPform.Y.H = 0;
        QPform.Y.lb = 0;
        QPform.Y.ub = inf;
    end
end

function [QPform, dX_dt,SSi] = loadElectricGenerator(Gen)
% this function loads the parameters for an electric generator generators
[QPform,dX_dt,SSi] = loadCHPGenerator(Gen);

function [QPform,dX_dt,SSi] = loadChiller(Gen)
% this function loads the parameters for a chiller generators
% it is very similar to the way Electric Generators are loaded
[QPform,dX_dt,SSi] = loadCHPGenerator(Gen);

function [QPform,dX_dt,SSi] = loadCHPGenerator(Gen)
% this function loads the parameters for a combined heat and power
% generator, regular electric generator, or chiller
global Plant
UB = Gen.Size;
if isfield(Gen.Output,'Cooling')&& Gen.Output.Cooling(end)>0
    LB = Gen.VariableStruct.Startup.Cooling(1); %chiller
    costTerms = GenCosts(Gen,LB,UB,'Cooling');
    if Plant.optimoptions.sequential
        Cratio = Gen.Output.Cooling(end);
        QPform.output.C = 1;
        QPform.output.E = -1/Cratio;
        costTerms.P = UB;
        costTerms.I = UB;
        costTerms.Intercept(1) =0;
        costTerms.Convex(1) = 0;
    else
        QPform.output.C = 1;
    end
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
if license('test','Control_System_Toolbox')
    [dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,Hratio);
else dX_dt = Gen.Size/10; SSi = [];
end
dX_dt = dX_dt/Plant.optimoptions.scaletime;    
QPform.Ramp.A = [-1, 1; 1, -1;]; %-output1+output2=ramp up %output1-output2=-rampdown
QPform.Ramp.b = [dX_dt;dX_dt];
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
QPform.Ramp.A = [-1, 1; 1, -1;];
QPform.Ramp.b = [dX_dt; dX_dt];

QPform.X.H = 0;
QPform.X.f = costLinear;
QPform.X.lb = 0;
QPform.X.ub = Gen.Size;

QPform.states(1,2) = {'X'};
QPform.X.H(2) = 0;%fit B
QPform.X.f(2) = costLinear;%fit B
QPform.X.lb(2) = LB; %fit B
QPform.X.ub(2) = Gen.Size;%fit B

function [QPform,dX_dt,SSi] = loadElectricStorage(Gen)
%this function loads all the parameters needed for electric storage
[QPform,dX_dt,SSi] = Storage(Gen);

function [QPform,dX_dt,SSi] = loadThermalStorage(Gen)
%loads all types of thermal storage
[QPform,dX_dt,SSi] = Storage(Gen);

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
        QPform.output.C = []; 
    elseif strcmp(Gen.VariableStruct.EnStoreType, 'HotTES')
        QPform.output.H = [];
        Stor.ChargeEff = 1; %set to ideal 
        Stor.DischEff = 1;
    elseif strcmp(Gen.VariableStruct.EnStoreType, 'HVAC')
        if strcmp(Gen.VariableStruct.EnStoreType, 'Heat only')
            QPform.output.H = [];
        elseif strcmp(Gen.VariableStruct.EnStoreType, 'AC only')
            QPform.output.C = [];
        else
            QPform.output.C = []; 
            QPform.output.H = [];
        end
    end
    
else %electric battery
    QPform.output.E = []; 
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
Outs = fieldnames(QPform.output);
dX_dt = Stor.PeakDisch;% storage discharge constraint in kW
a = (1/Stor.ChargeEff - Stor.DischEff);
if length(Outs)==2 %HVAC system
    QPform.states = {'X'; 'Y'}; %state of charge, charging power, no buffers
    QPform.link.ineq = [a -1]; % (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
    QPform.link.bineq = 0;
    
    QPform.X.lb = 0;
    QPform.X.ub = Stor.UsableSize;
    QPform.X.H = 0;
    QPform.X.f = 0;
    QPform.Ramp.A = [-1,1;1,-1]; % SOC2-SOC1<peakcharge; SOC1-SOC2<peakdischarge
    QPform.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];%since the stored energy is applied directly to the demand site, the max discharge rate is infinite.

    QPform.Y.lb = 0;
    QPform.Y.ub = Stor.PeakCharge;%the limit on how much charging power can be delivered is handled by the generators' limits, so put inf here to prevent redundancy
    QPform.Y.H = 0;%the cost of the charging power is handled by the generators 
    QPform.Y.f = 0;
elseif a==1 %ideal storage, ignore charging state
    QPform.states = {'X';'Z';'W'};%SOC(t+1), charging power, upper buffer, lower buffer
    QPform.link.ineq = [-1 0 -1; 1 -1 0];%-SOC(t)-lowerbuffer<-0.2  and %SOC-upperbuffer<0.8
    QPform.link.bineq = [0;0;]; %% note: the magnitude of the buffer is set later in findBuffer
    QPform.Stor.ChargeEff = 1; %eliminate chargin inefficiencies in CHP plants
    QPform.Stor.DischEff = 1;
    
    QPform.X.lb = 0;
    QPform.X.ub = Stor.UsableSize;
    QPform.X.H = 0; %no cost directly associated with storage
    QPform.X.f = 0;
    QPform.Ramp.A = [-1, 1; 1, -1];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
    QPform.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];

    QPform.Z.lb = 0;
    QPform.Z.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    QPform.Z.H = 0;
    QPform.Z.f = 0;

    QPform.W.lb = 0;
    QPform.W.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    QPform.W.H = 0;
    QPform.W.f = 0;
else % include charging state
    QPform.states = {'X';'Y';'Z';'W'};%SOC(t+1), charging power, upper buffer, lower buffer
    QPform.link.ineq = [a, -1, 0, 0];% (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
    QPform.link.ineq = [QPform.link.ineq; -1 0 0 -1];%-SOC(t)-lowerbuffer<-0.2
    QPform.link.ineq = [QPform.link.ineq; 1 0 -1 0];%SOC-upperbuffer<0.8
    QPform.link.bineq = [0;0;0;];
    
    QPform.X.lb = 0;
    QPform.X.ub = Stor.UsableSize;
    QPform.X.H = 0; %no cost directly associated with storage
    QPform.X.f = 0;
    QPform.Ramp.A = [-1, 1; 1, -1];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
    QPform.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];

    QPform.Y.lb = 0; %the lower bound of the charging state is 0 (only shows up as load when charging
    QPform.Y.ub = Stor.PeakCharge; %the upper bound of the charging power is limited by the generators themselves. 
    QPform.Y.H = 0;
    QPform.Y.f = 0;

    QPform.Z.lb = 0;
    QPform.Z.ub = 0;
    QPform.Z.H = 0;
    QPform.Z.f = 0;

    QPform.W.lb = 0;
    QPform.W.ub = 0;
    QPform.W.H = 0;
    QPform.W.f = 0;
end

function [QPform,dX_dt,SSi] = loadSolar(Gen)
%PV solar
QPform = [];
QPform.output.E = [];%there are no states or outputs for solar because renewable outputs are handled on the demand side
QPform.states = [];
dX_dt = inf;
SSi = [];

function [QPform, dX_dt,SSi] = loadHydroStorage(Gen)
% this function loads the parameters for a hydroelectric plant.
global Plant
SSi =[];
QPform.Stor.Size = Gen.Size*Plant.optimoptions.scaletime;
QPform.Stor.SelfDischarge  = 0; %needs to be evaporative losses
QPform.Stor.UsableSize  = QPform.Stor.Size*((Gen.VariableStruct.MaxHead-Gen.VariableStruct.MinHead)/Gen.VariableStruct.MaxHead);
Eff = Gen.VariableStruct.MaxGenCapacity/(Gen.VariableStruct.MaxGenFlow*Gen.VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW
QPform.output.E = Eff*Gen.VariableStruct.MaxHead*84.674;%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 87.674 kJ/ (1000ft^3*ft)
QPform.output.W = 1;
QPform.states = {'X';'Z';'W'};

QPform.X.H = 0;
QPform.X.f = 0;
QPform.X.lb = 0;
QPform.X.ub = Gen.Size;
QPform.Ramp.b = [Gen.VariableStruct.RampDown; Gen.VariableStruct.RampUp;];

%%buffer states
QPform.link.ineq = [-1 0  -1];%-SOC(t)-lowerbuffer<-0.2
QPform.link.ineq = [QPform.link.ineq; 1 -1 0];%SOC-upperbuffer<0.8
QPform.link.bineq = [0;0;];

QPform.Z.lb = 0;
QPform.Z.ub = 0;
QPform.Z.H = 0;
QPform.Z.f = 0;

QPform.W.lb = 0;
QPform.W.ub = 0;
QPform.W.H = 0;
QPform.W.f = 0;

dX_dt = Gen.VariableStruct.RampUp;