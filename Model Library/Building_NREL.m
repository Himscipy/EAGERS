function Plant = Building_NREL
global SimSettings
SimSettings.NominalTemperature = 22;
%% Components
Plant.Components.Building.type = 'BuildingNREL'; %fresh air 
Plant.Components.Building.name = 'Building';
Plant.Components.Building.CapacitanceRoom = 5e7;
Plant.Components.Building.CapacitanceWall = 5e7;
Plant.Components.Building.Resistance = 0.001;
Plant.Components.Building.ResistanceWall = 0.001;
Plant.Components.Building.Area = 313.42;
Plant.Components.Building.Occupancy = 5.38195521e-2;%people/m^2
Plant.Components.Building.OccupancySchedule = [0 0 0 0 0 0 .1 .2 .95 .95 .95 .95 .5 .95 .95 .95 .95 .3 .1 .1 .1 .1 .05 .05];
Plant.Components.Building.PlugLoad = 8.07007; %W/m^2
Plant.Components.Building.PlugSchedule = [.4 .4 .4 .4 .4 .4 .4 .4 .9 .9 .9 .9 .8 .9 .9 .9 .9 .5 .4 .4 .4 .4 .4 .4];
Plant.Components.Building.LightingLoad = 9.6875; % W/m^2
Plant.Components.Building.LightingSchedule = [.05 .05 .05 .05 .05 .1 .09238 .27714 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .4619 .27714 .27714 .18476 .18476 .09238 .04619];
Plant.Components.Building.TagInf = {'Temperature';'WallTemperature';};
Plant.Components.Building.connections = {'HVAC.massflow';'HVAC.temperature';'AmbientTemperature';};

%% Controls (note: controls can have specification that depends upon a initialized variable of a component)
Plant.Controls.HVAC.type = 'HVAC_NREL';
Plant.Controls.HVAC.name = 'HVAC';
Plant.Controls.HVAC.Target = {22.2}; %target dry bulb temperature 
Plant.Controls.HVAC.COP = 3.5;
Plant.Controls.HVAC.Entrainment = 40; % percent of fresh air
Plant.Controls.HVAC.Gain = 1e-2;%gain for T control
Plant.Controls.HVAC.PropGain = 10;
Plant.Controls.HVAC.TagInf = {'FlowRate';'Temperature';};
Plant.Controls.HVAC.connections = {'';'Building.Temperature';'Building.Mode';'AmbientTemperature';};

Plant.Scope = {'Building.Temperature';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'Building.WallTemperature';'HVAC.Temperature';'HVAC.FlowRate';}];
end%Ends function BuildingNREL