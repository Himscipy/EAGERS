function block = InitializeBuildingNREL(varargin)
% a simple building model
% Three (3) inlets: mass flow and temperature from HVAC system,  ambient temperature
% Three (3) outlets: Temperature, Wall Temperature, Mode (heating/cooling)
% Two (2) states: Room Temperature, Wall Temperature
global Tags
block = varargin{1};
if length(varargin)==1 % first initialization
    block.Scale = [22.2 23]; %temperature of room and wall
    block.IC = ones(length(block.Scale),1);
    
    block.Cp = 1e3; % J/kg*K
    
    %%
    block.InletPorts = {'massflow','temperature','Tamb'};
    block.massflow.IC = 4;
    block.temperature.IC = 12.8;
    block.Tamb.IC = 25;
    
    block.OutletPorts = {'Temperature';'WallTemperature';'Mode';};
    block.Temperature.IC = 22.2;
    block.WallTemperature.IC = 0.0085;
    block.Mode.IC = 'cooling';
    
    block.P_Difference = {};
    
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).WallTemperature = block.WallTemperature.IC;
    Tags.(block.name).Mode = 'cooling';
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    t = 0;
    h = mod(t/3600,24);
    T = block.Scale;
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)
    error = 1;
    while abs(error)>1e-3
        T(1) = 0.33*T(1) + 0.67*((T(2)/block.Resistance + IntGains + Inlet.massflow*block.Cp*Inlet.temperature)/(Inlet.massflow*block.Cp + 1/block.Resistance)); %temperature which  balances energy equation
        T(2) = 0.33*T(2) + 0.67*((T(1)/block.Resistance + Inlet.Tamb/block.ResistanceWall)/(1/block.Resistance + 1/block.ResistanceWall)); %temperature which  balances energy equation
        error = (block.Scale(1) - T(1));
        block.Scale = T; %temperature state is in celcius
    end
    
    block.Temperature.IC = block.Scale(1);
    block.WallTemperature.IC = block.Scale(2);
    if ((Inlet.Tamb - T(1))/block.Resistance + IntGains)>0
        block.Mode.IC = 'cooling';
    else
        block.Mode.IC = 'heating';
    end
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).WallTemperature = block.WallTemperature.IC;
    Tags.(block.name).Mode = block.Mode.IC;
end