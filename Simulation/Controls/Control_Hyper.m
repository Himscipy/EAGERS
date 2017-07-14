function Out = Control_Hyper(varargin)
% Controls for recoperated micro-turbine with hot bypass, cold bypass, and bleed valves
% Four (4) inlets: Two targets (Power, turbine exhaust temperature setpoint), Turbine exit temperature (TET), RPM, Combustor equivelance
% Five (5) outlets: Two corresponding to the targets, Generator Power, fuel flow rate, hot byass valve postion, cold bypass valve position, and compressor bleed valve position
% Five (5) states: dgenerator power, molar fuel flow, cold bypass, hot bypass, and bleed valve
% Like all controllers, if there are n targets, the first n inlets set 
% those targets, while the first n outlets correspond to the controllers 
% ability to hit that target
% The model can thus be linearized around any of these targets
global Tags 
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Load control';'Fuel Control';'Cold Bypass';'Hot Bypass';'Bleed'};
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

    block.OutletPorts(end+1:end+5) = {'GenPower','FuelFlow','ColdBypass','HotBypass','BleedValve'};
    block.GenPower.IC = block.NominalPower/block.GenEfficiency; 
    block.FuelFlow.IC = block.Target1.IC/(block.EstimatedEfficiency*(HeatingValue(block.Fuel))); %molar flow rate of fuel inlet
    block.ColdBypass.IC = 0;
    block.HotBypass.IC = 0;
    block.BleedValve.IC = 0;

    block.Scale = [block.NominalPower/block.GenEfficiency; block.FuelFlow.IC; 1; 1; 1];
    block.IC = [1; 1; block.ColdBypass.IC;block.HotBypass.IC; block.BleedValve.IC;];% inital condition 
    block.UpperBound = [inf,inf,1,1,1];
    block.LowerBound = [0,0,0,0,0];
    block.InitializeError = 1;
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    block.RPM.IC = Inlet.RPM;
    block.TET.IC = Inlet.TET;
    
    
    SteadyPower = ComponentProperty(strcat(block.SteadyPower,'.IC'));
    PowerError =(Inlet.Target1 - SteadyPower*block.GenEfficiency)/block.Target1.IC;
    block.InitializeError = abs(PowerError);
    
    if abs(PowerError)<0.06
        a = .3;
    else a = .15;
    end
    block.GenPower.IC = Inlet.Target1/block.GenEfficiency;
    block.FuelFlow.IC = block.FuelFlow.IC + a*PowerError*block.Scale(3); %adjust fuel flow to control power;
    block.Measured1.IC = block.GenPower.IC*block.GenEfficiency;
    block.Measured2.IC = Inlet.TET;
    block.ColdBypass.IC = block.IC(3);
    block.HotBypass.IC = block.IC(4);
    block.BleedValve.IC = block.IC(5);
    
    block.IC(1) = block.GenPower.IC/(block.Target1.IC/block.GenEfficiency);%adjust initial Power
    block.IC(2) = block.FuelFlow.IC/block.Scale(3); %adjust fuel flow to control TET
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    
    P_Gain = block.PropGain.*block.Scale;

    PowerError = Inlet.Target1-Y(2);
    GenPower =(P_Gain(2)*PowerError + Y(2));
    FuelError =(Inlet.Target1 - GenPower*block.GenEfficiency)/Inlet.Target1;

    if strcmp(string1,'Outlet')
        Out.Measured1.IC = GenPower*block.GenEfficiency;
        Out.Measured2.IC = Inlet.TET;
        Out.GenPower = GenPower;
        Out.FuelFlow = (P_Gain(2)*FuelError + Y(2));
        Out.ColdBypass = Y(3);
        Out.HotBypass = Y(4);
        Out.BleedValve = Y(5);
        Tags.(block.name).TET = Inlet.TET;
        Tags.(block.name).GenPower = Out.GenPower*block.GenEfficiency;
        Tags.(block.name).FuelFlow = Out.FuelFlow;
        Tags.(block.name).Efficiency = GenPower/(Out.FuelFlow*(HeatingValue(block.Fuel))); %molar flow rate of fuel inlet
    elseif strcmp(string1, 'dY')
        dY = Y*0;
        dY(1) = PowerError*block.IntGain(2);
        dY(2) = FuelError*block.IntGain(3);
        dY(3) = 0;
        Tags.(block.name).dY1 = PowerError;
        Tags.(block.name).dY2 = FuelError;
        Out = dY;
    end
end
end%Ends function Control_Hyper