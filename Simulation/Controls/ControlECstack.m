function Out = ControlECstack(varargin)
% Controls for electrolyzer stack only, control air flow to anode, inlet temperature and steam flow rate
% Four (4) inlets: Two Targets (temperature set point and power), Measured Temperature, Voltage
% Seven (7) outlets: Two measurements of targets, OxidentTemp, oxidant flow rate, steam temperature, steam flow rate, Current
% One  (1) state: If there is no oxidant flow, then control current to maintain temperature
% Like all controllers, if there are n targets, the first n inlets set 
% those targets, while the first n outlets correspond to the controllers 
% ability to hit that target
% The model can thus be linearized around any of these targets
global Tags
F=96485.339; % %Faraday's constant in Coulomb/mole
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Oxidant Flow Rate';'Fuel Cell Current';};
    block.TargetDescription = {'Operating Temperature';'Net Power'};
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
    block.InletPorts(end+1:end+2) = {'Temperature','Voltage'};
    block.Temperature.IC = ComponentProperty(strcat(block.connections{3},'.IC'));
    block.Voltage.IC = ComponentProperty(strcat(block.connections{4},'.IC'));
    
    block.deltaTStack = ComponentProperty(block.deltaTStack);
    block.Cells = ComponentProperty(block.Cells);
    block.Steam = ComponentProperty(block.Steam);
    block.Utilization = ComponentProperty(block.Utilization);
    block.SteamTemperature = ComponentProperty(block.SteamTemperature);
    
    OxFlow = NetFlow(ComponentProperty(block.OxidantFlow)); 
    Current = block.Target2.IC*1000/(block.Voltage.IC*block.Cells);
    SteamFlow = (block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    [h,~] = enthalpy(block.Target1.IC,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    block.Vbalance = 1./(2*F)*h_rxn3; %voltage that balances heat

    block.OutletPorts = {'OxidantTemp','OxidantFlow','SteamTemp','SteamFlow','Current'};
    block.OxidantTemp.IC = block.Target1.IC;
    block.OxidantFlow.IC = OxFlow; 
    block.SteamTemp.IC = block.SteamTemperature;
    block.SteamFlow.IC = SteamFlow;
    block.Current.IC = Current;
    
    block.P_Difference = {};

    block.IC = 1; % inital condition
        
    if OxFlow>0
        block.HasFlow = true;
        block.Scale = [OxFlow;];
        block.UpperBound = inf;
        block.LowerBound = 0;
    else 
        block.HasFlow = false;
        block.Scale = [Current;];
        block.UpperBound = inf;
        block.LowerBound = -inf;
    end
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    PowerError = (Inlet.Target2-abs(block.Current.IC)*Inlet.Voltage*block.Cells/1000)/Inlet.Target2;
    block.Current.IC = block.Current.IC*(1 + PowerError);
    block.SteamFlow.IC = (block.Cells*abs(block.Current.IC)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    block.SteamTemp.IC = block.SteamTemp.IC; %currently not manipulated during initialization
    
    averageT = mean(Inlet.Temperature);
    
    block.Measured1.IC = averageT;
    block.Measured2.IC = block.Current.IC*Inlet.Voltage*block.Cells/1000; %Power in kW
    
    if block.HasFlow
        Q_cathode = block.OxidantFlow.IC*40*block.deltaTStack;
        
        if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(block.Current.IC)/1000) - Q_cathode)>0
            block.OxidantTemp.IC = Inlet.Target1-100; %cooling stack
            TavgError = (averageT -Inlet.Target1)/block.deltaTStack; % too hot = increase flow
        else
            block.OxidantTemp.IC = Inlet.Target1+100;%heating stack
            TavgError = (Inlet.Target1-averageT)/block.deltaTStack; %too hot = reduce flow
        end
        block.InitializeError = abs(TavgError); 

        newFlow = block.OxidantFlow.IC*(1+.5*TavgError);
        block.OxidantFlow.IC = newFlow;
        block.Scale = [block.OxidantFlow.IC];
        block.IC = [1-TavgError*block.PropGain(1)]; %Oxidant flow rate
    else
        block.InitializeError = 0;
        block.Scale = [block.Current.IC];
        block.IC = [1-PowerError*block.PropGain(2)]; % inital condition with no anode flow in: net current
    end
    Out = block;
else %running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    averageT = mean(Inlet.Temperature);
    % averageT = (Inlet.Hot - Inlet.Cold + block.dT_cath_PEN);

    if block.HasFlow
        Current = -Inlet.Target2*1000/(Inlet.Voltage*block.Cells);
    else
        Current = Y(1);
        VoltError = block.Vbalance - Inlet.Voltage;
    end
    SteamFlow = block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000;

    if block.HasFlow
        Q_cathode = Y(1)*40*block.deltaTStack;
        if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(Current)/1000) - Q_cathode)>0
            Out.OxidantTemp = Inlet.Target1-100; %cooling stack
            TavgError = (averageT -Inlet.Target1)/block.deltaTStack; % too hot = increase flow
        else
            Out.OxidantTemp = Inlet.Target1+100;%heating stack
            TavgError = (Inlet.Target1-averageT)/block.deltaTStack; %too hot = reduce flow
        end
        Out.OxidantFlow = max(0,Y(1)+(TavgError*block.PropGain(1))*block.Scale(1));
    else
        Out.OxidantTemp = Inlet.Target1;
        Out.OxidantFlow = 0;
    end

    if strcmp(string1,'Outlet')
        Out.Measured1.IC = averageT;
        Out.Measured2.IC = Current*Inlet.Voltage*block.Cells/1000; %Power in kW
        Out.SteamTemp = block.SteamTemp.IC;%currently no control of fuel inlet temp
        Out.SteamFlow = SteamFlow;
        Out.Current = Current;
        Tags.(block.name).OxidantTemp = Out.OxidantTemp;
        Tags.(block.name).OxidantFlow = Out.OxidantFlow;
        Tags.(block.name).Tsteam = block.SteamTemp.IC;%currently no control of steam inlet temp
        Tags.(block.name).SteamFlow = SteamFlow;
        Tags.(block.name).Current = Current;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        if block.HasFlow
            if Y(1)>0 || TavgError>0 %anti-windup so flow does not go negative
                dY(1) = block.Gain(1)*TavgError;
            end
        else
            dY(1) = block.Gain(2)*VoltError;
        end
        Out = dY;
    end
end
end%Ends function ControlECstack