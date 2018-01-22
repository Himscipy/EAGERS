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
if any(strcmp(Gen.Source, {'NG';'Diesel';}))%if it is a fuel utility
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
        QPform.X.lb = 0;
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
n = 5; %%need to select n based on how good the fit is.
[QPform,dX_dt,SSi] = loadPiecewise(Gen,n);
end%Ends function loadElectricGenerator

function [QPform, dX_dt,SSi] = loadCHPGenerator(Gen)
% this function loads the parameters for an electric generator generators
n = 5; %%need to select n based on how good the fit is.
[QPform,dX_dt,SSi] = loadPiecewise(Gen,n);
end%Ends function loadCHPGenerator

function [QPform,dX_dt,SSi] = loadChiller(Gen)
% global Plant
n = 5; %%need to select n based on how good the fit is.
[QPform,dX_dt,SSi] = loadPiecewise(Gen,n);
end%Ends function loadChiller

function [QPform,dX_dt,SSi] = loadHeater(Gen)
n = 5; %%need to select n based on how good the fit is.
[QPform,dX_dt,SSi] = loadPiecewise(Gen,n);
end%Ends function loadHeater

function [QPform,dX_dt,SSi] = loadPiecewise(Gen,n)
% this function loads the parameters for a combined heat and power
% generator, regular electric generator, or chiller
% n is number of segments
% order is either 1 or 2, when order is 2 it solves for the quadratic 
% coefficients of C = c_0 + a_1*x_1 + a_2*x_2 + ... a_n*X_n + b_1*x_1^2 + b_2*x_2^2 + ... b_n*X_n^2
% subject to b_i >0, and a_i> a_(i-1) + b_(i-1)*(x_i)_max
% if order is 1 it solves for linear coefficients 
% of C = c_0 + a_1*x_1 + a_2*x_2 + ... a_n*X_n
%subject to  a_i> a_(i-1)
global Plant
UB = Gen.Size;
order = 2;
if isfield(Gen.Output,'Cooling')&& Gen.Output.Cooling(end)>0
    order = 1;
    LB = Gen.VariableStruct.Startup.Cooling(end); %chiller
    QPform.output.C = 1;
    out = 'Cooling';
elseif isfield(Gen.Output,'Electricity')&& Gen.Output.Electricity(end)>0
    LB = Gen.VariableStruct.Startup.Electricity(end); %electric or CHP generator 
    QPform.output.E = 1;
    out = 'Electricity';
elseif isfield(Gen.Output,'Heat')&& Gen.Output.Heat(end)>0
    LB = Gen.VariableStruct.Startup.Heat(end);
    QPform.output.H = 1;
    out = 'Heat';
end
% in Fit A: c_0 = 0
capacity = Gen.Output.Capacity;
efficiency = Gen.Output.(out);
operationRange = find(capacity>=LB/UB);
[P,I] = sort(capacity(operationRange));
eff = efficiency(operationRange);
Y = P./eff(I); %cost of generator in terms of input at outputs P
Y(isnan(Y)) = 0;

if n ==2 && (isfield(Gen.Output,'Electricity') || isfield(Gen.Output,'Cooling'))
    n2 = 10;
    costA = zeros(n2,1);
    costB = zeros(n2,1);
    for i = 1:1:n2
        ub = [LB/UB + (UB-LB)/UB*(i-1)/n2,1];
        [costA(i),~,~,~] = fitfcn(P,Y,ub,n,order,0);
        [costB(i),~,~,~] = fitfcn(P,Y,ub,n,order,1);
    end
    [~,I_a] = min(costA);
    ub = [LB/UB + (UB-LB)/UB*(I_a-1)/n2,1];
    [~,A,~,~] = fitfcn(P,Y,ub,n,order,0);
    [~,I_b] = min(costB);
    ub = [LB/UB + (UB-LB)/UB*(I_b-1)/n2,1];
    [~,B,c_0,x_max] = fitfcn(P,Y,ub,n,order,1);
else
    ub = LB/UB + (UB-LB)/UB*(1:n)/n;
    [~,A,~,~] = fitfcn(P,Y,ub,n,order,0);
    [~,B,c_0,x_max] = fitfcn(P,Y,ub,n,order,1);
end

%% scale by size of gen (later can normalize by changing QPform.Output so that all states are 0-1)
[cap,I] = sort(capacity);
x_max = x_max*UB;
c_0 = c_0*UB;
cap = cap*UB;
if order == 2
    A(n+1:2*n) = A(n+1:2*n)/UB;
    B(n+1:2*n) = B(n+1:2*n)/UB;
end

if isfield(Gen.Output,'Electricity')&& Gen.Output.Electricity(end)>0 && isfield(Gen.Output,'Heat')&& Gen.Output.Heat(end)>0
    Hratio_0 = Gen.Output.Heat(I)./Gen.Output.Electricity(I);
    Hratio_0(isinf(Hratio_0)) = 0;
    QPform.output.H(1,1) = interp1(cap,Hratio_0,min(x_max(1),mean([LB,x_max(1)])));
    for i = 2:1:n
        QPform.output.H(i,1) = min(QPform.output.H(i-1,1),mean(interp1(cap,Hratio_0,linspace(max(LB,sum(x_max(1:i-1))),min(UB,sum(x_max(1:i))),10))));
    end
    Hratio_vector = zeros(n,1);
    Xcum = zeros(n+1,1);
    Xcum(1) = LB;
    Xcum(2) = x_max(1);
    for i = 2:1:n
        Xcum(i+1) = Xcum(i)+x_max(i);
    end
    for i = 1:1:n
        Hratio_vector(i) = (interp1(cap,cap.*Hratio_0,Xcum(i+1)) - interp1(cap,cap.*Hratio_0,Xcum(i)))/(Xcum(i+1)-Xcum(i));
    end
    H_0 = interp1(cap,cap.*Hratio_0,LB) - LB*Hratio_vector(1);
    
    Hratio_midway = interp1(Xcum,[0;Hratio_vector],(LB+UB)/2);
    QPform.output.H(1,2) = max(Hratio_vector(1),Hratio_midway);
    for i = 2:1:n
        if i<n/2
            Hratio_vector(i) = max(Hratio_vector(i),Hratio_midway);%once past midway allow Hratio to drop below midway value
        end
        QPform.output.H(i,2) =  min(QPform.output.H(i-1,2),Hratio_vector(i));
    end
    QPform.constDemand.H = -H_0;
    x_max_A = x_max;
    x_max_B = x_max;
else %remove segments that are equal (same slope & quadratic) to previous segment
    keep = true(1,order*n);
    x_max_A = x_max(1);
    for i = 2:1:n
        if A(i)-A(i-1)<1e-5 && (order == 1 || abs(A(n+i)-A(n+i-1)) < 1e-10)
            keep(i) = false;
            if order ==2
                keep(n+i) = false;
            end
            x_max_A(end) = x_max_A(end) + x_max(i);
        else
            x_max_A(end+1) = x_max(i);
        end
    end
    x_max_A = x_max_A*UB/sum(x_max_A);
    x_max_A(end) = UB-sum(x_max_A(1:end-1)); %make sure that there is no rounding error, and the last x_max_A value actually meets the UB
    A = A(keep);

    keep = true(1,order*n);
    x_max_B = x_max(1);
    for i = 2:1:n
        if B(i)-B(i-1)<1e-5 && (order == 1 || abs(B(n+i)-B(n+i-1)) < 1e-10)
            keep(i) = false;
            if order == 2
                keep(n+i) = false;
            end
            x_max_B(end) = x_max_B(end) + x_max(i);
        else
            x_max_B(end+1) = x_max(i);
        end
    end
    x_max_B = x_max_B*UB/sum(x_max_B);
    x_max_B(end) = UB-sum(x_max_B(1:end-1));%make sure there is no rounding error
    B = B(keep);
end
n_A = length(A)/order;
n_B = length(B)/order;

%put into form that buildMatrices will see
letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';};
QPform.states = cell(n,2);
QPform.states(1:n_A,1) = letters(1:n_A);
QPform.states(1:n_B,2) = letters(1:n_B);
if strcmp(out,'Cooling')
    if strcmp(Gen.Source,'Electricity')
        source = 'E';
    else
        source = 'H';
    end
    QPform.output.(source)(1:n_A,1) = -A(1:n_A);
    QPform.output.(source)(1:n_B,2) = -B(1:n_B);
    QPform.constDemand.(source) = c_0;
    A = A*0;
    B = B*0;
elseif strcmp(out,'Electricity') || strcmp(out,'Heat')
    QPform.constCost = c_0;
end

for i = 1:1:n
    if i<=n_A
        if order == 2
            QPform.(letters{i}).H = A(n_A+i);
        else QPform.(letters{i}).H = 0;
        end
        QPform.(letters{i}).f = A(i);
        QPform.(letters{i}).lb = 0;
        QPform.(letters{i}).ub = x_max_A(i);
    end
    if i<=n_B
        if order == 2
            QPform.(letters{i}).H(2) = B(n_B+i);
        else QPform.(letters{i}).H(2) = 0;
        end
        QPform.(letters{i}).f(2) = B(i);
        QPform.(letters{i}).ub(2) = x_max_B(i);
        if i ==1
            QPform.(letters{i}).lb(2) = LB;
            QPform.(letters{i}).ub(2) = max(LB,x_max_B(i));
        else
            QPform.(letters{i}).lb(2) = 0;
        end
    end
end

% % %% plot
% figure
% plot(P*UB,Y*UB,'b')
% hold on
% X2 = linspace(0,x_max_A(1),10);
% if order == 1
%     Y2 = linspace(0,x_max_A(1),10)*A(1);
% else Y2 = linspace(0,x_max_A(1),10)*A(1) + linspace(0,x_max_A(1),10).^2*A(n_A+1);
% end
% for i = 2:1:n_A
%     X2 = [X2, X2(end)+linspace(0,x_max_A(i),10)];
%     if order == 1
%         Y2 = [Y2,Y2(end)+linspace(0,x_max_A(i),10)*A(i)];
%     else Y2 = [Y2,Y2(end)+linspace(0,x_max_A(i),10)*A(i) + linspace(0,x_max_A(i),10).^2*A(n_A+i)];
%     end
% end
% plot(X2,Y2,'g')
% 
% X3 = linspace(0,x_max_B(1),10);
% if order == 1
%     Y3 = c_0+linspace(0,x_max_B(1),10)*B(1);
% else Y3 = c_0+linspace(0,x_max_B(1),10)*B(1) + linspace(0,x_max_B(1),10).^2*B(n_B+1);
% end
% for i = 2:1:n_B
%     X3 = [X3, X3(end)+linspace(0,x_max_B(i),10)];
%     if order == 1
%         Y3 = [Y3,Y3(end)+linspace(0,x_max_B(i),10)*B(i)];
%     else Y3 = [Y3,Y3(end)+linspace(0,x_max_B(i),10)*B(i) + linspace(0,x_max_B(i),10).^2*B(n_B+i)];
%     end
% end
% plot(X3,Y3,'r')

%% ramp rates
if isfield(Gen.Output,'Electricity')&& Gen.Output.Electricity(end)>0 && isfield(Gen.Output,'Heat')&& Gen.Output.Heat(end)>0
    H_0 = QPform.output.H(1,end)*LB;
else H_0 = [];
end
if license('test','Control_Toolbox')
    [dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,H_0);
else dX_dt = Gen.Size/1; SSi = []; %assume 4 hours from off to peak
end
dX_dt = dX_dt/Plant.optimoptions.scaletime;  
QPform.Ramp.b = [dX_dt;dX_dt]; %-output1+output2=ramp up %output1-output2=-rampdown
end%Ends function loadPiecewise

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
a = (1/Stor.ChargeEff - Stor.DischEff);
QPform.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
if a~=0 %not ideal storage, add charging state
    QPform.states = {'X'; 'Y'}; %state of charge, charging power, no buffers
    QPform.link.ineq = [a -1]; % (1/nc - nd)*(SOC(t) -SOC(t-1)) - charging <0 ------ Charging is the 1/inefficiency + 1
    QPform.link.bineq = 0;
    QPform.Y.lb = 0;
    QPform.Y.ub = Stor.PeakCharge;%the limit on how much charging power can be delivered is handled by the generators' limits, so put inf here to prevent redundancy
    QPform.Y.H = 0;%the cost of the charging power is handled by the generators 
    QPform.Y.f = 0;
end    
if isfield(Gen.VariableStruct,'Buffer') && Gen.VariableStruct.Buffer ~= 0 %buffer states
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
QPform.Stor.DischEff = 1; %100% effecient
%Put this at 1/2 full capacity until the MaxHead-MinHead/MaxHead is fixed
%Currently, with the previous solution, some plants lost a vast majority of usable space
% QPform.Stor.UsableSize  = QPform.Stor.Size**((Gen.VariableStruct.MaxHead-Gen.VariableStruct.MinHead)/Gen.VariableStruct.MaxHead);
month = datevec(Plant.Data.Timestamp(1));
month = month(2); %month starting in
QPform.Stor.UsableSize = QPform.Stor.Size*Gen.VariableStruct.StartingPoint(1,month);
QPform.Stor.Power2Flow = 1/(Eff*Gen.VariableStruct.MaxHead*84.674);%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 84.674 kJ/ (1000ft^3*ft)
QPform.output.E = [1;0;];
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
    QPform.output.E = [1; 0; 0;];
end
if isfield(Gen.VariableStruct,'Buffer') && Gen.VariableStruct.Buffer ~= 0 %buffer states
    if Gen.VariableStruct.MaxSpillFlow>0
        QPform.states = {'X';'Y';'S';'U';'L'};
        QPform.link.eq = [QPform.Stor.Power2Flow, 0, 1, 0, 0];
        QPform.link.ineq = [0 -1 0 0  -1; 0 1 0 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        QPform.output.E = [1; 0; 0; 0; 0;];
    else
        QPform.states = {'X';'Y';'U';'L'};
        QPform.link.eq = [QPform.Stor.Power2Flow, 0, 0, 0];
        QPform.link.ineq = [0 -1 0  -1; 0 1 -1 0;];%-SOC(t)-lowerbuffer<-0.2, %SOC-upperbuffer<0.8
        QPform.output.E = [1; 0; 0; 0;];
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

function [cost,A,c_0,x_max] = fitfcn(P,Y,ub,n,order,intercept)
x_max = zeros(1,n);
x_max(1) = ub(1);
x_max(2:n) = ub(2:n) - ub(1:n-1);
X = zeros(length(P),n);
for i = 1:1:length(P)
    X(i,1) = min(P(i),ub(1));
    for j = 2:1:n
        X(i,j) = max(0,min(x_max(j),P(i)-sum(X(i,1:j-1))));
    end
end
QP.H = ones(order*n);
QP.f = ones(order*n,1);
for i = 1:1:order*n
    if i<=n
        a = X(:,i);
    else
        a = X(:,i-n).^2;
    end
    for j = 1:1:order*n
        if j<=n
            b = X(:,j);
        else
            b = X(:,j-n).^2;
        end
        QP.H(i,j) = 2*sum(a.*b);
    end
    QP.f(i) = -2*sum(a.*Y);
end

QP.A = zeros(order*n-1,order*n);
QP.b = zeros(order*n-1,1);
for i = 1:1:n-1
    QP.A(i,i) = 1;
    QP.A(i,i+1) = -1;
    if order==2
        QP.A(i,n+i) = 2*x_max(i);
    end
end
if order==2
    for i = n:1:2*n-1
        QP.A(i,i+1) = -1;
    end
    QP.Aeq = [x_max x_max.^2];
else
    QP.Aeq = x_max;
end
QP.beq = Y(end);
QP.lb = [];
QP.ub =[];
QP.solver = 'quadprog';
if intercept
    a = zeros(1,order*n);
    for i = 1:1:order*n
        if i<=n
            a(i) = 2*sum(X(:,i));
        else
            a(i) = 2*sum(X(:,i-n).^2);
        end
    end
    QP.H = [2*length(P), a; [a',QP.H];];
    QP.f = [-2*sum(Y);QP.f];
    QP.A = [zeros(order*n-1,1),QP.A];
    QP.Aeq = [1 QP.Aeq];
end
A = callQPsolver(QP);
cost = 0.5*A'*QP.H*A + A'*QP.f;
if order == 2
    for i = n+1:2*n
        if A(i) <1e-6
            A(i) = 0;
        end
    end
end
if intercept
    c_0 = A(1);
    A = A(2:end);
else c_0 = 0;
end
end%Ends fitfcn