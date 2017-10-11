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
    OxFlow = NetFlow(ComponentProperty(block.OxidantFlow)); 
    if OxFlow>0
        block.PIdescription = {'Oxidant Flow Rate';'Oxidant Temperature';};
    else block.PIdescription = {'Current'};
    end
    block.TargetDescription = {'Operating Temperature';'Hot/Cold Temperature Differential';'Net Power'};
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
        block.(strcat('Target',num2str(i))).Saturation =  [-inf,inf];
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+2) = {'Hot','Voltage'};
    block.Hot.IC = ComponentProperty(strcat(block.connections{4},'.IC'));
    block.Hot.Saturation =  [0,inf];
    block.Voltage.IC = ComponentProperty(strcat(block.connections{5},'.IC'));
    block.Voltage.Saturation =  [-inf,inf];
    
    block.Cells = ComponentProperty(block.Cells);
    block.Steam = ComponentProperty(block.Steam);
    block.Utilization = ComponentProperty(block.Utilization);
    block.SteamTemperature = ComponentProperty(block.SteamTemperature);
    
    Current = block.Target3.IC*1000/(block.Voltage.IC*block.Cells);
    SteamFlow = (block.Cells*abs(Current)/(2000*F*block.Utilization*block.Steam.H2O));
    [h,~] = enthalpy(block.Target1.IC,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    block.Vbalance = h_rxn3./(2*F); %voltage that balances heat

    block.OutletPorts = {'OxidantTemp','OxidantFlow','SteamTemp','SteamFlow','Current'};
    block.OxidantTemp.IC = block.Target1.IC;
    block.OxidantFlow.IC = OxFlow; 
    block.SteamTemp.IC = block.SteamTemperature;
    block.SteamFlow.IC = SteamFlow;
    block.Current.IC = Current;
    
    block.P_Difference = {};

    if OxFlow>0
        block.IC = [1, 1]; % inital condition
        block.HasFlow = true;
        block.Scale = [OxFlow;block.OxidantTemp.IC];
        block.UpperBound = inf;
        block.LowerBound = 0;
    else
        block.IC = 1; % inital condition
        block.HasFlow = false;
        block.Scale = Current;
        block.UpperBound = inf;
        block.LowerBound = -inf;
    end
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    PowerError = (Inlet.Target3-block.Current.IC*Inlet.Voltage*block.Cells/1000)/Inlet.Target3;
    block.Current.IC = block.Current.IC*(1 + PowerError);
    block.SteamFlow.IC = (block.Cells*abs(block.Current.IC)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    block.SteamTemp.IC = block.SteamTemp.IC; %currently not manipulated during initialization
    
    PEN_Temperature = mean(ComponentProperty('EC1.T.Elec'));
    averageT = (mean(Inlet.Hot)+block.OxidantTemp.IC)/2; %average temperature of cathode 
    block.dT_cath_PEN = PEN_Temperature - averageT; %temperature differencce between cathode and PEN 
    
    
    averageT = (block.OxidantTemp.IC+mean(Inlet.Hot))/2+block.dT_cath_PEN;
    
    block.Measured1.IC = averageT;
    block.Measured2.IC = block.Current.IC*Inlet.Voltage*block.Cells/1000; %Power in kW
    
    if block.HasFlow
        deltaT = (mean(Inlet.Hot)-block.OxidantTemp.IC);
        dTerror =(deltaT-Inlet.Target2)/Inlet.Target2;

        TavgError = (Inlet.Target1-averageT)/Inlet.Target2; 
        block.InitializeError = max(abs(TavgError),abs(dTerror));

        a = .5;
        block.OxidantFlow.IC = block.OxidantFlow.IC*(1+a*dTerror); 
        block.OxidantTemp.IC = block.OxidantTemp.IC + a*(TavgError + .75*dTerror)*Inlet.Target2;
        
        block.Scale = [block.OxidantFlow.IC;block.OxidantTemp.IC];
        block.IC = [1-dTerror*block.PropGain(1);1-TavgError*block.PropGain(2)]; %Oxidant flow rate
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
    Inlet = checkSaturation(Inlet,block);
    P_Gain = block.PropGain.*block.Scale;
    if block.HasFlow
        Current = Inlet.Target3*1000/(Inlet.Voltage*block.Cells);
    else
        Current = Y(1);
        VoltError = block.Vbalance - Inlet.Voltage;
    end
    SteamFlow = block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000;

    if block.HasFlow
        averageT = (Y(2)+mean(Inlet.Hot))/2+block.dT_cath_PEN;
        TavgError = (Inlet.Target1-averageT)/Inlet.Target2; %too hot = reduce flow
        deltaT = abs(mean(Inlet.Hot)-(Y(2)+TavgError*P_Gain(2)));
        dTerror =(deltaT-Inlet.Target2)/Inlet.Target2;
        OxFlow = Y(1)+dTerror*P_Gain(1);
        OxTemp = Y(2)+TavgError*P_Gain(2);
        if OxFlow < block.OxidantFlow.IC*block.MinOxidantFlowPerc/100
            OxFlow = block.OxidantFlow.IC*block.MinOxidantFlowPerc/100;
            FixedFlow = true;
        else
            FixedFlow = false;
        end
    else
        OxTemp = Inlet.Target1;
        OxFlow = 0;
    end

    if strcmp(string1,'Outlet')
        Out.Measured1.IC = averageT;
        Out.Measured2.IC = Current*Inlet.Voltage*block.Cells/1000; %Power in kW
        Out.OxidantFlow = OxFlow;
        Out.OxidantTemp = OxTemp;
        Out.SteamTemp = block.SteamTemp.IC;%currently no control of fuel inlet temp
        Out.SteamFlow = SteamFlow;
        Out.Current = Current;
        Tags.(block.name).OxidantTemp = Out.OxidantTemp;
        Tags.(block.name).OxidantFlow = Out.OxidantFlow;
        Tags.(block.name).Tsteam = block.SteamTemp.IC;%currently no control of steam inlet temp
        Tags.(block.name).SteamFlow = SteamFlow;
        Tags.(block.name).Current = Current;
        Tags.(block.name).deltaT = deltaT;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        Gain = block.Gain.*block.Scale;
        if block.HasFlow
            if FixedFlow %avoid reducing air flow as power goes to zero
                dY(1) = Gain(1)*(OxFlow - Y(1));
            else
                dY(1) = Gain(1)*dTerror;
            end
            dY(2) = Gain(2)*TavgError;
        else
            dY(1) = Gain(1)*VoltError;
        end
        Out = dY;
    end
end
end%Ends function ControlECstack