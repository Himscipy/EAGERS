function Plant = Building3
%multi-zone building
%% Components
Plant.Components.Building.type = 'ZonalBuilding'; %fresh air 
Plant.Components.Building.name = 'Building';
Plant.Components.Building.Capacitance = 5e7; %capacitance of each zone
Plant.Components.Building.Resistance = 0.001;
Plant.Components.Building.Area = 313.42;  %area of each zone
Plant.Components.Building.Volume = 940;
Plant.Components.Building.maxHumidity = 0.012;
%occupance, plug load and lighting can be specified 1 for the entire
%building, or seperately for each zone. Each row is a different zone
Plant.Components.Building.Occupancy = 5.38195521e-2;%people/m^2
Plant.Components.Building.OccupancySchedule = [0 0 0 0 0 0 .1 .2 .95 .95 .95 .95 .5 .95 .95 .95 .95 .3 .1 .1 .1 .1 .05 .05];
Plant.Components.Building.PlugLoad = 8.07007; %W/m^2
Plant.Components.Building.PlugSchedule = [.4 .4 .4 .4 .4 .4 .4 .4 .9 .9 .9 .9 .8 .9 .9 .9 .9 .5 .4 .4 .4 .4 .4 .4];
Plant.Components.Building.LightingLoad = 9.6875; % W/m^2
Plant.Components.Building.LightingSchedule = [.05 .05 .05 .05 .05 .1 .09238 .27714 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .4619 .27714 .27714 .18476 .18476 .09238 .04619];
Plant.Components.Building.TagInf = {'Temperature';'Humidity';'Cooling';'Heating';};
Plant.Components.Building.connections = {'AmbientTemperature';'AmbientHumidity';'HVAC.AirFlow';'HVAC.Tset';'HVAC.DPset';'HVAC.Damper';'HVAC.ColdWater';'HVAC.HotWater';};

%% Controls (note: controls can have specification that depends upon a initialized variable of a component)
Plant.Controls.HVAC.type = 'HVACproportional';
Plant.Controls.HVAC.name = 'HVAC';
Plant.Controls.HVAC.zones = 1;
Plant.Controls.HVAC.Target = [22.2; 11;]; %target dry bulb and dew point temperatures 
Plant.Controls.HVAC.Tolerance = 0.5; % deviation in temperature that triggers chiller
Plant.Controls.HVAC.Entrainment = 40; % percent of fresh air
Plant.Controls.HVAC.ColdAirSetpoint = 12.8; %minimum temperature of cooling air supplied by HVAC
Plant.Controls.HVAC.HotAirSetpoint = 40; %temperature of heating air supplied by HVAC when furnace is on
Plant.Controls.HVAC.minFlow = 4;
Plant.Controls.HVAC.maxFlow = 16;
Plant.Controls.HVAC.COP = 3.5;
Plant.Controls.HVAC.Gain = [1e-2; 0; 1e-2; 1e-2;];% integral gain for fan speed, damper position, coold water flow, hot water flow
Plant.Controls.HVAC.PropGain = [10; 0; 0; 0;];
Plant.Controls.HVAC.TagInf = {'FlowRate';'Temperature';};
Plant.Controls.HVAC.connections = {'Building.Temperature';'AmbientTemperature';'Building.Qcool';'Building.Qheat';'Building.Mode';};

Plant.Scope = {'Building.Temperature';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'Building.Humidity';'HVAC.Temperature';'HVAC.FlowRate';'HVAC.Cooling';'HVAC.Heating';}];