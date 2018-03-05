function locateBuildings
%identify what subNet node of electrical heating and cooling the building is on, and if there are chillers/heaters connected, disengage the internal HVAC
global Plant
build = [];
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
    buildLocation = cell(nB,1);
    for i = 1:1:nB
        buildLocation(i,1) = {Plant.Building(i).Location};
    end
    networkNames = {'Electrical';'DistrictHeat';'DistrictCool';};
    for net = 1:1:length(networkNames)
        for n = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
            node_m = Plant.subNet.(networkNames{net}).nodes{n};
            for k = 1:1:length(node_m)
                for i = 1:1:nB
                    if any(strcmp(buildLocation{i,1},node_m))
                        build(end+1) = i;
                        if strcmp(networkNames{net},'Electrical')
                            Plant.Building(i).QPform.Electrical.subnetNode = n;
                            Plant.Building(i).QPform.Location = Plant.subNet.Electrical.Location(n);
                        elseif strcmp(networkNames{net},'DistrictHeat')
                            Plant.Building(i).QPform.DistrictHeat.subnetNode = n;
                            if ~isempty(Plant.subNet.DistrictHeat.Equipment{n})
                                equip = Plant.subNet.DistrictHeat.Equipment{n};
                                for j = 1:1:length(equip)
                                    if ismember(Plant.Generator(equip(j)).Type,{'Heater';'CHP Generator';})
                                        Plant.Building(i).QPform.Heating = true;
                                    end
                                end
                            end
                            if ~isempty(Plant.subNet.DistrictHeat.connections{n}) %connected to heaters at a different node
                                Plant.Building(i).QPform.Heating = true;
                            end
                            if Plant.Building(i).QPform.Heating
                                Plant.Building(i).QPform.H2E = Plant.Building(i).VariableStruct.FanPower/(1.025*(50-15))/1.225; %Flow = Heating/(Cp_Air*(Tsupply(t)-Tmix)); Power = FanPower*Flow/air density
                            end
                        elseif strcmp(networkNames{net},'DistrictCool')
                            Plant.Building(i).QPform.DistrictCool.subnetNode = n;    
                            if ~isempty(Plant.subNet.DistrictCool.Equipment{n})
                                equip = Plant.subNet.DistrictCool.Equipment{n};
                                for j = 1:1:length(equip)
                                    if ismember(Plant.Generator(equip(j)).Type,{'Chiller';})
                                        Plant.Building(i).QPform.Cooling = true;
                                    end
                                end
                            end
                            if ~isempty(Plant.subNet.DistrictHeat.connections{n}) %connected to heaters at a different node
                                Plant.Building(i).QPform.Cooling = true;
                            end
                            if Plant.Building(i).QPform.Cooling
                                Plant.Building(i).QPform.C2E = Plant.Building(i).VariableStruct.FanPower/(1.025*(25-12))/1.225; %Flow = Cooling/(Cp_Air*(Tmix-Tsupply(t))); Power = FanPower*Flow/air density
                            end
                        end
                    end
                end
            end
            Plant.subNet.(networkNames{net}).Buildings{n} = build;
        end
    end
end
end%Ends function locateBuildings