function Plant = Chiller1
global SimSettings
SimSettings.NominalCooling = 2250*3.517;
%% Components
Plant.Components.CoolingTowerPump.type = 'WaterPump';
Plant.Components.CoolingTowerPump.name = 'CoolingTowerPump';
Plant.Components.CoolingTowerPump.HeadLoss = 10; % m of head loss that pump must supply
Plant.Components.CoolingTowerPump.Efficiency = .75; 
Plant.Components.CoolingTowerPump.TagInf = {'Flow_GPS';};
Plant.Components.CoolingTowerPump.connections = {'Controller.CoolingTowerPump';'Chiller.CTreturnTemp';};

Plant.Components.ColdWaterPump.type = 'WaterPump';
Plant.Components.ColdWaterPump.name = 'ColdWaterPump';
Plant.Components.ColdWaterPump.HeadLoss = 20; % m of head loss that pump must supply
Plant.Components.ColdWaterPump.Efficiency = .75; 
Plant.Components.ColdWaterPump.TagInf = {'Flow_GPS';};
Plant.Components.ColdWaterPump.connections = {'Controller.ColdWaterPump';11}; %should be connected to building CEreturn temperature

Plant.Components.CoolingTower.type = 'CoolingTower';
Plant.Components.CoolingTower.name = 'CoolingTower';
Plant.Components.CoolingTower.hconv = 30;
Plant.Components.CoolingTower.Cooling_Tons = SimSettings.NominalCooling/3.517*1.2;
Plant.Components.CoolingTower.HeaterDeltaT = 2;
Plant.Components.CoolingTower.HeadLoss = .3; % m of head loss that fan must supply
Plant.Components.CoolingTower.FanEfficiency = .75;
Plant.Components.CoolingTower.Mass = 300;
Plant.Components.CoolingTower.CP = 0.5; %specific heat of solid material (kJ/kg*K)
Plant.Components.CoolingTower.TagInf = {'Temperature';'WetBulbTemperature'};
Plant.Components.CoolingTower.connections = {'Controller.CoolingTowerFan';'CoolingTowerPump.Flow';'AmbientTemperature';'AmbientHumidity';};

Plant.Components.Chiller.type = 'CentrifugalChiller'; 
Plant.Components.Chiller.name = 'Chiller';
Plant.Components.Chiller.Refrigerant = 'R123';
Plant.Components.Chiller.Refrig_Tons = SimSettings.NominalCooling/3.517; %Refrigerationin tons (3.517kW)
Plant.Components.Chiller.HeaterDeltaT = 2; %aproximate deltaT across heat exchangers
Plant.Components.Chiller.EstimatedCOP = 6;
Plant.Components.Chiller.hconv_ev = 30; %convective heat transfer in evaporator
Plant.Components.Chiller.hconv_cd = 25; %convective heat transfer in condensor
Plant.Components.Chiller.Volume_ev = 1; %refrigerant volume (m^3)
Plant.Components.Chiller.Volume_cd = 1;  %refrigerant volume (m^3)
Plant.Components.Chiller.CP_ev = 0.5; %specific heat of solid material (kJ/kg*K)
Plant.Components.Chiller.CP_cd = 0.5; %specific heat of solid material (kJ/kg*K)
Plant.Components.Chiller.Mass_ev = 300;
Plant.Components.Chiller.Mass_cd = 300;
Plant.Components.Chiller.compdesP = 115/35;
Plant.Components.Chiller.RPMdesign = 1e4;
Plant.Components.Chiller.FlowDesign = Plant.Components.Chiller.Refrig_Tons*.025;%Mass flow of refrigerant in kg/s use a factor of roughly 0.0244kg/s per ton of chilling for R123
Plant.Components.Chiller.ShaftDensity = 1000; %kg/m^3
Plant.Components.Chiller.ShaftLength = 0.2; %m 
Plant.Components.Chiller.ShaftRadius = 0.2; %m
Plant.Components.Chiller.PeakEfficiency = .90;
Plant.Components.Chiller.MotorEfficiency = .90;
Plant.Components.Chiller.Map = 'CentrifugalChiller'; % Loads a saved compressor map
Plant.Components.Chiller.TagInf = {'RPM';'EvaporatorPressure';'CondensorPressure';'ColdWaterTemp';'Cooling';};
Plant.Components.Chiller.connections = {'Controller.Chiller1';'ColdWaterPump.Flow';'CoolingTower.Return';};

%% Controls (note: controls can have specification that depends upon a initialized variable of a component)
Plant.Controls.Controller.type = 'SingleChiller';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = [4; 22;]; %target chilled water temperatures (chiller/cooling tower)
Plant.Controls.Controller.NominalCoolingCapacity = SimSettings.NominalCooling;
Plant.Controls.Controller.EstimatedCOP = 5;
Plant.Controls.Controller.Gain = [1e-4; 1e-4; 1e-4; 1e-4;];%integral gain for power control
Plant.Controls.Controller.PropGain = [.5 .5 .5 .5];%proportional gain for power control
Plant.Controls.Controller.TagInf = {'Power';'COP'};
Plant.Controls.Controller.connections = {'CoolingDemand';11;'Chiller.ColdWaterTemp';'Tags.ColdWaterPump.Flow_GPS';'Chiller.CTreturnTemp';'CoolingTower.Temperature';'Tags.CoolingTowerPump.Flow_GPS';};

Plant.Scope = {'Chiller.Cooling';'Chiller.RPM';'CoolingTower.WetBulbTemperature'}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'CoolingTower.Temperature';'Controller.NetPower';'Controller.COP';}];

SimSettings.CoolingTime = [0 24];
SimSettings.CoolingDemand = [1 1]*Plant.Components.Chiller.Refrig_Tons*3.517;
SimSettings.TemperatureTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
SimSettings.TemperatureData = 24*ones(1,25);
SimSettings.HumidityTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
SimSettings.HumidityData = 0.01*ones(1,25);