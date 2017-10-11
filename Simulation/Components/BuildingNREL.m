function Out = BuildingNREL(varargin)
% a simple building model
% Three (3) inlets: mass flow and temperature from HVAC system,  ambient temperature
% Three (3) outlets: Temperature, Wall Temperature, Mode (heating/cooling)
% Two (2) states: Room Temperature, Wall Temperature
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.Scale = [22.2 23]; %temperature of room and wall
    block.IC = ones(length(block.Scale),1);
    block.UpperBound = [inf,inf];
    block.LowerBound = [0,0];
    
    block.Cp = 1e3; % J/kg*K
    
    %%
    block.InletPorts = {'massflow','temperature','Tamb'};
    block.massflow.IC = 4;
    block.massFlow.Saturation = [0,inf];
    block.temperature.IC = 12.8;
    block.temperature.Saturation = [10,30];
    block.Tamb.IC = 25;
    block.Tamb.Saturation = [-50,50];
    
    block.OutletPorts = {'Temperature';'WallTemperature';'Mode';};
    block.Temperature.IC = 22.2;
    block.WallTemperature.IC = 0.0085;
    block.Mode.IC = 'cooling';
    
    block.P_Difference = {};
    
%     Tags.(block.name).Temperature = block.Temperature.IC;
%     Tags.(block.name).WallTemperature = block.WallTemperature.IC;
%     Tags.(block.name).Mode = 'cooling';
   Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    
    Inlet = checkSaturation(Inlet,block);
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
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    
    Inlet = checkSaturation(Inlet,block);
    h = mod(t/3600,24);
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)

    if strcmp(string1,'Outlet')
        Out.Temperature = Y(1);
        Out.WallTemperature = Y(2);
        if ((Inlet.Tamb - Y(1))/block.Resistance + IntGains)>0
            Out.Mode = 'cooling';
        else Out.Mode = 'heating';
        end
        Tags.(block.name).Temperature = Y(1);
        Tags.(block.name).WallTemperature = Y(2);
        Tags.(block.name).Mode = Out.Mode;
    elseif strcmp(string1,'dY')
        dY = 0*Y;
        dY(1) =(((Y(2) - Y(1))/block.Resistance) + IntGains + Inlet.massflow*block.Cp*(Inlet.temperature - Y(1)))/block.CapacitanceRoom;%Change in temperature
        dY(2) = (((Y(1) - Y(2))/block.Resistance) + ((Inlet.Tamb - Y(2))/block.ResistanceWall))/block.CapacitanceWall;%change in wall temperature
        Out = dY;
    end
end%Ends function BuildingNREL