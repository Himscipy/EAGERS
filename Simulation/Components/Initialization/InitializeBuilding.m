function block = InitializeBuilding(varargin)
% a simple building model
% Two (2) inlets: air flow from HVAC system,  ambient temperature
% Three (3) outlets: Temperature, Humidity, Mode (heating/cooling)
% Two (2) state: Temperature, Humidity
global Tags
block = varargin{1};
if length(varargin)==1 % first initialization
    block.Scale = [22.2 .0085]; %temperature and humidity
    block.IC = ones(length(block.Scale),1);
    %%
    block.InletPorts = {'HVAC','Tamb'};
    block.HVAC.IC = makeAir(12.8,50,4,'rel');
    block.Tamb.IC = 25;
    
    block.OutletPorts = {'Temperature';'Humidity';'Mode';};
    block.Temperature.IC = 22.2;
    block.Humidity.IC = 0.0085;
    block.Mode.IC = 'cooling';
    
    block.P_Difference = {};
    
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
    Tags.(block.name).Mode = 'cooling';
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    t = 0;
    h = mod(t/3600,24);
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    Hin = (MassFlow(Inlet.HVAC) - 18*Inlet.HVAC.H2O)*enthalpyAir(Inlet.HVAC); %mass flow of dry air* enthalpy of dry air
    OutFlow = Inlet.HVAC;
    OutFlow.H2O = OutFlow.H2O + Occupancy*3e-6;%3x10^-6 kg/s is the average rate of water vapor released per person due to respiration  

    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)
    error = 1;
    while abs(error)>1e-3
        Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O)*enthalpyAir(OutFlow);
        Herror = 1/(1000*block.Resistance)*(Inlet.Tamb + block.Resistance*(IntGains + Hin*1000 - Hout*1000) + 273.15 - OutFlow.T); %temperature which  balances energy equation
        OutFlow.T = OutFlow.T + Herror/(1.5*MassFlow(OutFlow));
        error = (block.Scale(1) - (OutFlow.T-273.15));
        block.Scale(1) = OutFlow.T-273.15; %temperature state is in celcius
    end
    
    block.Scale(2) = OutFlow.H2O*18/(MassFlow(OutFlow)-OutFlow.H2O*18);%absolute humidity
    block.Temperature.IC = block.Scale(1);
    block.Humidity.IC = block.Scale(2);
    if ((Inlet.Tamb - (OutFlow.T-273.15))/block.Resistance + IntGains)>0
        block.Mode.IC = 'cooling';
    else
        block.Mode.IC = 'heating';
    end
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
    Tags.(block.name).Mode = block.Mode.IC;
end