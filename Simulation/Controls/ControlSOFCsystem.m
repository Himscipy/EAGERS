function Out = ControlSOFCsystem(varargin)
% Controls for SOFC system, control blower power, heater bypass, anode recirculation and fuel flow rate
% Tries to control the FC inlet temperature using heater bypass (too hot) and fuel utilization (too cold)
% Tries to control the FC exit temperature using the air flow indirectly through blower power
% Measured voltage determines current and fuel flow to meet power and utilization
% Seven (7) Inlets: Four Targets (temperature set point, temperature differential, steam to carbon ration, and power), T FC exit, T FC inlet, Voltage
% Nine (9) outlets: Four associated with the targets, Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Two (2) states: Heater bypass, blower power
% May need to add state for current back in to avoid fuel starvation during step changes
global Tags F
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Heater Bypass';'Blower Power';'Fuel Cell Current';};
    block.TargetDescription = {'Operating Temperature';'Hot/Cold Temperature Differential';'Steam to Carbon Ratio';'Net Power'};
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target  = Target;
    
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+3) = {'Hot','Cold','Voltage'};
    block.Hot.IC = ComponentProperty(strcat(block.connections{5},'.IC'));
    block.Cold.IC = ComponentProperty(strcat(block.connections{6},'.IC'));
    block.Voltage.IC = ComponentProperty(strcat(block.connections{7},'.IC'));
    
    block.OutletPorts(end+1:end+5) = {'HeaterBypass','Blower','Current','AnodeRecirc','FuelFlow'};
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);
    block.Utilization = ComponentProperty(block.Utilization);

    Bypass = ComponentProperty(block.HeaterBypass);
    Blower = ComponentProperty(block.Blower);
    BlowerMass = ComponentProperty('Blower.FlowDesign');
    CathodeMass = NetFlow(ComponentProperty('FC1.Flow2Out.IC'))*28.84;
    ComponentProperty('Blower.FlowDesign',CathodeMass);
    Blower = Blower*CathodeMass/BlowerMass; %scale blower power
    Current = block.Target4.IC*1000/(block.Voltage.IC*block.Cells);
    Recirculation = ComponentProperty(block.AnodeRecirc);
    FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);

    block.HeaterBypass.IC = Bypass;
    block.Blower.IC = Blower; 
    block.Current.IC =  Current;
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.IC = FuelFlow;
    
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target3.IC;
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            block.WGSeffective = effectiveWGS(block.Fuel,FuelFlow,0.7,Steam2Carbon,block.Cells*Current/(2*F*1000),Recirculation);
        end
    end
    
    block.P_Difference = {};
    
    block.IC = [Bypass;1;]; % inital condition
    block.Scale = [1;Blower;];
    Out = block;
elseif length(varargin)==2 %second initialization
    block = varargin{1};
    Inlet = varargin{2};
    PEN_Temperature = mean(ComponentProperty('FC1.T.Elec'));
    BlowerPower = ComponentProperty('Blower.NominalPower');
    StackPower = Inlet.Target4 + BlowerPower;
%     PowerError = (StackPower-block.Current.IC*Inlet.Voltage*block.Cells/1000)/StackPower;
%     block.Current.IC = block.Current.IC*(1 + PowerError);
    
    block.Current.IC = StackPower*1000/(Inlet.Voltage*block.Cells);
    averageT = (mean(Inlet.Hot)+mean(Inlet.Cold))/2; %average temperature of cathode
    block.dT_cath_PEN = PEN_Temperature - averageT; %temperature differencce between cathode and PEN
    TinletError = (Inlet.Target1 - (mean(Inlet.Cold) + .5*Inlet.Target2 + block.dT_cath_PEN))/Inlet.Target2;
    
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-(Inlet.Target2))/Inlet.Target2;

    Steam2Carbon = Inlet.Target3;
    
    %% adjust blower mass flow to get deltaT correct
    blowerMassFlow = ComponentProperty('Blower.FlowDesign');
    ComponentProperty('Blower.FlowDesign',blowerMassFlow*(1+.5*dTerror));
    block.Blower.IC = BlowerPower*(1+.5*dTerror);
%     blower = block.Blower.IC
    %%if too cold, add more fuel, it will burn and more heat will be recovered
%     Utilization = block.Utilization*(1-0.04*TavgError)
%     block.Utilization = Utilization;
    
%     %% adjust heat exchanger effectiveness to control inlet temperature
%     HXtarget = ComponentProperty('HX1.Effectiveness');
%     HXtarget = min(0.96,HXtarget*(1 + .05*TavgError));
%     ComponentProperty('HX1.Target',HXtarget);
%     ComponentProperty('HX1.sizemethod','Effectiveness');
    
    %% adjust heat exchanger cold exit temperature to control FC inlet temperature
    HXtarget = ComponentProperty('HX1.Target');
    HXtarget = HXtarget + TinletError*Inlet.Target2;
    ComponentProperty('HX1.Target',HXtarget);

    %%adjust exit temperature to avoid energy feedback during iterations
    Terror = (Inlet.Target1 - (averageT+block.dT_cath_PEN));
    CathodeOut = ComponentProperty('FC1.Flow2Out.IC');
    CathodeOut.T = CathodeOut.T + Terror;
    ComponentProperty('FC1.Flow2Out.IC',CathodeOut);
    AnodeOut = ComponentProperty('FC1.Flow1Out.IC');
    AnodeOut.T = AnodeOut.T + Terror;
    ComponentProperty('FC1.Flow1Out.IC',AnodeOut);

    FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    block.FuelFlow.IC = FuelFlow;
    
    if block.AnodeRecirc.IC > 0
        block.AnodeRecirc.IC = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*block.Current.IC/(2*F*1000),block.AnodeRecirc.IC);
    end
    block.Measured1.IC = averageT;
    block.Measured2.IC = deltaT;
    block.Measured3.IC = Steam2Carbon;
    block.Measured4.IC = StackPower - BlowerPower;
    
    block.InitializeError = max(abs(TinletError),abs(dTerror)); 
%     block.InitializeError = abs(TavgError);
    block.Scale = [1; block.Blower.IC;];
    block.IC = [block.HeaterBypass.IC-TinletError*block.PropGain(1); 1-dTerror*block.PropGain(2);];    
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    I_Gain = block.Gain.*block.Scale;
    TinletError = ((mean(Inlet.Cold) + .5*Inlet.Target2 + block.dT_cath_PEN) - Inlet.Target1)/Inlet.Target2; %target a fixed inlet temperature
    ToutletError = ((mean(Inlet.Hot) - .5*Inlet.Target2 + block.dT_cath_PEN) - Inlet.Target1)/Inlet.Target2; %target a fixed outlet temperature
    averageT = (mean(Inlet.Hot) + mean(Inlet.Cold))/2;
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    % dTerror =(deltaT-Inlet.Target2)/Inlet.Target2; %target a fixed temperature gradient
    % BlowerPower = Y(2)+(dTerror*block.PropGain(2))*block.Scale(2);

    BlowerPower = Y(2)+(ToutletError*block.PropGain(2))*block.Scale(2);
    StackPower = Inlet.Target4 + BlowerPower;
    % Current = Y(3)+block.PropGain(3)*(StackPower*1000/(Inlet.Voltage*block.Cells) - Y(3));
    % Power = Current*Inlet.Voltage*block.Cells/1000;
    % PowerError = (StackPower - Power)/StackPower;

    Current = StackPower*1000/(Inlet.Voltage*block.Cells);

    Bypass = min(1,max(0,Y(1)+TinletError*block.PropGain(1)));
    Utilization = min(block.Utilization,block.Utilization + .1*(Y(1)+TinletError*block.PropGain(1))); %lower utilization if additional pre-heating is necessary
    FuelFlow = block.Cells*Current/(2*F*Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;

    if block.AnodeRecirc.IC == 0
        Recirculation = 0;
    else
        Steam2Carbon = Inlet.Target3;
        Recirculation = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*Current/(2*F*1000),block.AnodeRecirc.IC);
    end

    if strcmp(string1,'Outlet')
        Out.Measured1 = averageT;
        Out.Measured2 = deltaT;
        Out.Measured3 = Steam2Carbon;
        Out.Measured4 = StackPower - BlowerPower;
        Out.HeaterBypass = Bypass;
        Out.Blower = BlowerPower;
        Out.Current = Current;
        Out.AnodeRecirc = Recirculation;
        Out.FuelFlow = FuelFlow;
        Tags.(block.name).Bypass = Bypass;
        Tags.(block.name).Blower = BlowerPower;
        Tags.(block.name).Recirculation = Recirculation;
        Tags.(block.name).FuelFlow = FuelFlow;
        Tags.(block.name).Current = Current;
        Tags.(block.name).Utilization = Utilization;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        if Y(1)+TinletError*block.PropGain(1) < 0 
            dY(1) = 0.02*I_Gain(1)*TinletError; %controlling utilization
        else
            dY(1) = I_Gain(1)*TinletError;
        end
        %% need to deterimine saturation on blower
        if Y(1) <0
            dY(2) = 0;
        else
            dY(2) = I_Gain(2)*ToutletError;
        end

    %     dY(3) = I_Gain(3)*PowerError;
        Out = dY;  
    end
end
end%Ends function ControlSOFCsystem