function Out = BuildingNREL(t,Y, Inlet,block,string1)
% a simple building model
% Three (3) inlets: mass flow and temperature from HVAC system,  ambient temperature
% Three (3) outlets: Temperature, Wall Temperature, Mode (heating/cooling)
% Two (2) states: Room Temperature, Wall Temperature
global Tags
if strcmp(string1,'Outlet')
    Out.Temperature = Y(1);
    Out.WallTemperature = Y(2);
    Out.Mode = Tags.(block.name).Mode;
    Tags.(block.name).Temperature = Y(1);
    Tags.(block.name).WallTemperature = Y(2);
elseif strcmp(string1,'dY')
    h = mod(t/3600,24);
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)

    dY = 0*Y;
    dY(1) =(((Y(2) - Y(1))/block.Resistance) + IntGains + Inlet.massflow*block.Cp*(Inlet.temperature - Y(1)))/block.CapacitanceRoom;%Change in temperature
    dY(2) = (((Y(1) - Y(2))/block.Resistance) + ((Inlet.Tamb - Y(2))/block.ResistanceWall))/block.CapacitanceWall;%change in wall temperature
    if ((Inlet.Tamb - Y(1))/block.Resistance + IntGains)>0
        Tags.(block.name).Mode = 'cooling';
    else Tags.(block.name).Mode = 'heating';
    end
    Out = dY;
end