function Out = ControlRecouperatedGasTurbine(varargin)
% Controls for recoperated micro-turbine
% Four (4) inlets: Two targets (Power, turbine exhaust temperature setpoint), Turbine exit temperature (TET), RPM, Combustor equivelance
% Four (4) outlets: Two corresponding to the targets, Generator Power, fuel flow rate
% Three (3) states: desired RPM, generator power, molar fuel flow
% Like all controllers, if there are n targets, the first n inlets set 
% those targets, while the first n outlets correspond to the controllers 
% ability to hit that target
% The model can thus be linearized around any of these targets
global Tags 
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Target RPM';'Power Pulled by the Generator';'Fuel Flow';};
    block.TargetDescription = {'Net Power';'Turbine Exit Temperature';};
    
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
    block.InletPorts(end+1:end+2) = {'TET','RPM'};
    block.TET.IC = ComponentProperty(strcat(block.connections{3},'.IC'));
    block.RPM.IC = ComponentProperty(strcat(block.connections{4},'.IC'));
    
    block.OutletPorts(end+1:end+2) = {'GenPower','FuelFlow'};
    block.GenPower.IC = block.Target1.IC/block.GenEfficiency; 
    block.FuelFlow.IC = block.Target1.IC/(block.EstimatedEfficiency*HeatingValue(block.Fuel)); %molar flow rate of fuel inlet
    
    block.InitializeError = 1;
    
    block.Scale = [block.RPMdesign; block.Target1.IC/block.GenEfficiency; block.FuelFlow.IC;];
    block.IC = [1; 1; 1;]; % inital condition 
    block.UpperBound = [2,inf,inf];
    block.LowerBound = [0,0,0];
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    block.RPM.IC = Inlet.RPM;
    block.TET.IC = Inlet.TET;
    
    block.GenPower.IC = Inlet.Target1/block.GenEfficiency;
%     TETerror = (Inlet.Target2-Inlet.TET)/400; %error in TET -- should change RPM
    SteadyPower = ComponentProperty(strcat(block.SteadyPower,'.IC'));
    PowerError =(Inlet.Target1 - SteadyPower*block.GenEfficiency)/block.Target1.IC;
    block.InitializeError = abs(PowerError);
    
    if abs(PowerError)<0.06
        a = .3;
    else a = .15;
    end
    block.FuelFlow.IC = block.FuelFlow.IC + a*PowerError*block.Scale(3); %adjust fuel flow to control power;
    block.Measured1.IC = block.GenPower.IC*block.GenEfficiency;
    block.Measured2.IC = Inlet.TET;
    
    block.IC(1) = Inlet.RPM/block.RPMdesign;%adjust initial RPM
    block.IC(2) = block.GenPower.IC/(block.Target1.IC/block.GenEfficiency);%adjust initial Power
    block.IC(3) = block.FuelFlow.IC/block.Scale(3); %adjust fuel flow to control TET
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    P_Gain = block.PropGain.*block.Scale;
    I_Gain = block.Gain.*block.Scale;

    TETerror = (Inlet.TET-Inlet.Target2)/100; 
    RPMfromTET =(P_Gain(1)*TETerror + Y(1));
    RPMerror = (Inlet.RPM - RPMfromTET)/block.RPMdesign;
    GenPower =(P_Gain(2)*RPMerror + Y(2));
    FuelError =(Inlet.Target1 - GenPower*block.GenEfficiency)/Inlet.Target1;
    PowerError = RPMerror;

    if strcmp(string1,'Outlet')
        Out.Measured1.IC = GenPower*block.GenEfficiency;
        Out.Measured2.IC = Inlet.TET;
        Out.GenPower = GenPower;
        Out.FuelFlow = Y(3) + P_Gain(3)*FuelError;
        Tags.(block.name).TET = Inlet.TET;
        Tags.(block.name).GenPower = Out.GenPower*block.GenEfficiency;
        Tags.(block.name).FuelFlow = Out.FuelFlow;
        Tags.(block.name).Efficiency = GenPower/(Out.FuelFlow*(HeatingValue(block.Fuel))); %molar flow rate of fuel inlet
    elseif strcmp(string1, 'dY')
        dY = Y*0;
        dY(1) = TETerror*I_Gain(1);
        dY(2) = PowerError*I_Gain(2);
        dY(3) = FuelError*I_Gain(3);
        Tags.(block.name).dY1 = TETerror;
        Tags.(block.name).dY2 = PowerError;
        Tags.(block.name).dY3 = FuelError;
        Out = dY;
    end
end
end%Ends function ControlRecouperatedGasTurbine