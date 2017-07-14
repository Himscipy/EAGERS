function Out = ControlFCstack(varargin)
% Controls for Fuel cell stack only. This controller modulates air flow, 
% inlet temperature, fuel flow rate, net current and anode recirculation
% to meet a target average stack temperature and a net power production
% For a air cooled fuel cell the oxidant temperature and flow are
% controlled with PI controllers, the current is calculated from the target 
% power and the measured voltage, and the fuel flow and recirculation are
% calculated from the current and the specified fuel composition
% For the oxyFC recirculation, fuel flow and current use PI controllers,
% and the oxidant flow is calculated from the current.
% Six (6) inlets: Four Targets (temperature set point, temperature differential, steam to carbon ration, and power), T FC exit, Voltage
% Ten (10) outlets: Four measurements of targets, OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate, current, anode recirculation
% One or Three (1-3) states: oxidant flow rate or recirculation, fuel flow and current
% Like all controllers, if there are n targets, the first n inlets set 
% those targets, while the first n outlets correspond to the controllers 
% ability to hit that target
% The model can thus be linearized around any of these targets
global Tags
F=96485.339; % %Faraday's constant in Coulomb/mole
if length(varargin)==1 %first initialization
    block = varargin{1};
    if isfield(block,'OxyFC')
        block.PIdescription = {'Anode Recirculation';'Fuel Flow';'Fuel Cell Current';};
    else
        block.PIdescription = {'Oxidant Flow Rate';};
    end
    block.TargetDescription = {'Operating Temperature';'Hot/Cold Temperature Differential';'Steam to Carbon Ratio';'Net Power'};
    %convert text or pull variables from other blocks
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(block.Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target = Target;
    
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+2) = {'Hot','Voltage'};
    block.Hot.IC = ComponentProperty(strcat(block.connections{5},'.IC'));
    block.Voltage.IC = ComponentProperty(strcat(block.connections{6},'.IC'));
    
    block.OutletPorts(end+1:end+5) = {'OxidantTemp','OxidantFlow','AnodeRecirc','FuelFlow','Current'};
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);
    Current = block.Target4.IC*1000/(block.Voltage.IC*block.Cells);
    Recirculation = ComponentProperty(block.AnodeRecirc);
    if isfield(block,'OxyFC')
        FuelFlow = NetFlow(ComponentProperty(block.FuelFlow));
        block.Oxidant = ComponentProperty(block.Oxidant);
        OxidantTemp = block.Target(1);
        block.OxidantUtilization = ComponentProperty(block.OxidantUtilization);
        OxidantFlow = Current*block.Cells/(4000*F*block.Oxidant.O2*block.OxidantUtilization);
    else
        block.Utilization = ComponentProperty(block.Utilization);
        FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        OxidantTemp = ComponentProperty(block.OxidantTemp);
        OxidantFlow = NetFlow(ComponentProperty(block.OxidantFlow));
    end
    block.OxidantTemp.IC = OxidantTemp;
    block.OxidantFlow.IC = OxidantFlow; 
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.IC = FuelFlow; 
    block.Current.IC = Current;
    
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target(3);
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            block.WGSeffective = effectiveWGS(block.Fuel,FuelFlow,0.7,Steam2Carbon,block.Cells*Current/(2*F*1000),Recirculation);
        end
    end

    if isfield(block,'OxyFC')
        block.IC = [Recirculation;1;1;]; % inital condition
        block.Scale = [1;FuelFlow;Current];
        block.UpperBound = [1,inf,inf];
        block.LowerBound = [0,0,0];
    else
        block.IC = [1]; % inital condition
        block.Scale =  [OxidantFlow;];
        block.UpperBound = inf;
        block.LowerBound = 0;
    end
    
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};

    Power = block.Current.IC*Inlet.Voltage*block.Cells/1000;
%     PEN_Temperature = mean(ComponentProperty('FC1.T.Elec'));
    averageT = (mean(Inlet.Hot)+block.OxidantTemp.IC)/2; %average temperature of cathode
%     block.dT_cath_PEN = PEN_Temperature - averageT; %temperature differencce between cathode and PEN    
    TavgError = (Inlet.Target1-averageT)/Inlet.Target2;
    deltaT = (mean(Inlet.Hot)-block.OxidantTemp.IC);
    dTerror =(deltaT-Inlet.Target2)/Inlet.Target2;
     
    Steam2Carbon = Inlet.Target3;
%     dTerror
%     TavgError
    if isfield(block,'OxyFC')
        FuelFlow = block.FuelFlow.IC*(1-.02*TavgError);
        O2flow = block.Cells*block.Current.IC/(4*F*1000)*32*3600*24/1000; %Ton/day
        Parasitic = (1.0101*O2flow^-.202)*(O2flow*1000/24); %parasitic in kW
        PowerError = (Inlet.Target4 + Parasitic - Power)/Inlet.Target4;
        block.Current.IC = block.Current.IC*(1 + PowerError);
%         block.Current.IC = (Inlet.Target4 + Parasitic)*1000/(Inlet.Voltage*block.Cells);
    else
        block.Current.IC = Inlet.Target4*1000/(Inlet.Voltage*block.Cells);
        FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        a = .5;
        block.OxidantFlow.IC = block.OxidantFlow.IC*(1+a*dTerror); 
        block.OxidantTemp.IC = block.OxidantTemp.IC + a*(TavgError + .75*dTerror)*block.Target(2);
    end
    
    block.FuelFlow.IC = FuelFlow;
    if block.AnodeRecirc.IC > 0
        block.AnodeRecirc.IC = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*block.Current.IC/(2*F*1000),block.AnodeRecirc.IC);
    end
    block.Measured1.IC = averageT;
    block.Measured2.IC = deltaT;
    block.Measured3.IC = Steam2Carbon;
    block.Measured4.IC = Power;
    
    if isfield(block,'OxyFC')
        block.InitializeError = 0;%max(abs(PowerError));
%         block.Scale = [1, block.FuelFlow.IC];
%         block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2)];%-VoltageError*block.PropGain(3)];
        block.Scale = [1, block.FuelFlow.IC, block.Current.IC];
        block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2),1];%-VoltageError*block.PropGain(3)];
    else
        block.InitializeError = max(abs(TavgError),abs(dTerror));
        block.Scale = [block.OxidantFlow.IC];
        block.IC = [1-dTerror*block.PropGain(1)];
    end
    
    Out = block;
else %running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    averageT = (block.OxidantTemp.IC+mean(Inlet.Hot))/2;
    TavgError = (Inlet.Target1-averageT)/Inlet.Target2;
    deltaT = (mean(Inlet.Hot)-block.OxidantTemp.IC);
    dTerror =(deltaT-Inlet.Target2)/Inlet.Target2;
    %     
    if isfield(block,'OxyFC')
        Current = Y(3)+ (TavgError*block.PropGain(3))*block.Scale(3); % current controller for managing temperature in closed-end cathode FC (proportional control only)
        O2flow = block.Cells*Current/(4*F*1000)*32*3600*24/1000; %Ton/day
        Parasitic = (1.0101*O2flow^-.202)*(O2flow*1000/24); %parasitic in kW
        Power = Current*Inlet.Voltage*block.Cells/1000;
        PowerError = (Inlet.Target4 + Parasitic - Power)/(Inlet.Target4 + Parasitic);

        FuelFlow = Y(2)+(PowerError*block.PropGain(2))*block.Scale(2);
        CH4_Util = Y(3)*block.Cells/(2000*F)/(FuelFlow*4*block.Fuel.CH4);%hydrogen used vs. ideal H2 available
        R.CH4 = Current/(8000*F*CH4_Util);
        a = 4352.2./Inlet.Target1 - 3.99;
        R.WGS = R.CH4*exp(a)*block.WGSeffective;% Water gas shift equilibrium constant
        R.H2 = Current/(2000*F); %# of electrochemical reactions
        [h,~] = enthalpy(averageT,{'H2','H2O','O2','CO','CO2','CH4'});
        h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
        h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
        h_rxn3 = h.H2O - h.H2 -.5*h.O2;

        FuelIn = block.Fuel;
        FuelIn.T = block.OxidantTemp.IC;
        Cp = SpecHeat(FuelIn);

        BalancedPower = R.H2*(-h_rxn3) - (R.CH4*h_rxn1+R.WGS*h_rxn2) - FuelFlow/block.Cells*Cp*Inlet.Target2;
        BalancedVoltage = BalancedPower*1000/Current;
        VoltageError = Inlet.Voltage - BalancedVoltage;
    else
        Current = Inlet.Target4*1000/(Inlet.Voltage*block.Cells);
        Power = Current*Inlet.Voltage*block.Cells/1000;
        FuelFlow = block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;
        if block.AnodeRecirc.IC == 0
            Recirculation = 0;
        else
            Steam2Carbon = Inlet.Target3;
            Recirculation = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*Current/(2*F*1000),block.AnodeRecirc.IC);
        end
    end
    if strcmp(string1,'Outlet')
        Out.Measured1 = averageT;
        Out.Measured2 = deltaT;
        Out.Measured3 = Steam2Carbon;
        Out.Measured4 = Power;
        if isfield(block,'OxyFC')
            Out.OxidantFlow = block.Cells*Current/(4000*F*block.Oxidant.O2); % kmol/s
            Out.AnodeRecirc = Y(1)+ dTerror*block.PropGain(1)*Inlet.Target2;
        else
            Out.OxidantFlow = Y(1)+(dTerror*block.PropGain(1))*block.Scale(1);
            Out.AnodeRecirc =  Recirculation;
        end
        Out.OxidantTemp = block.OxidantTemp.IC;
        Out.FuelFlow = FuelFlow;
        Out.Current = Current;
        Tags.(block.name).FuelFlow = Out.FuelFlow;
        Tags.(block.name).Current = Out.Current;
        Tags.(block.name).Recirculation = Out.AnodeRecirc;
        Tags.(block.name).OxidantTemp = Out.OxidantTemp;
        Tags.(block.name).OxidantFlow = Out.OxidantFlow;
        Tags.(block.name).Power = Power;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        if isfield(block,'OxyFC')
            dY(1) = block.Gain(1)*dTerror; %recirculation changes temperature gradient
            dY(2) = block.Gain(2)*PowerError;%fuel flow
            dY(3) = block.Gain(3)*VoltageError;%current
        else  
            dY(1) = block.Gain(1)*dTerror;
        end
        Out = dY;
    end
end
end %Ends function ControlFCstack