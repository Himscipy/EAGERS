function Building = locate_buildings(Building,subNet,Gen)
%identify what subNet node of electrical heating and cooling the building is on, and if there are chillers/heaters connected, disengage the internal HVAC
build = [];
nB = length(Building);
buildLocation = cell(nB,1);
for i = 1:1:nB
    buildLocation(i,1) = {Building(i).Location};
end
networkNames = {'Electrical';'DistrictHeat';'DistrictCool';};
for net = 1:1:length(networkNames)
    for n = 1:1:length(subNet.(networkNames{net}).nodes)
        node_m = subNet.(networkNames{net}).nodes{n};
        for k = 1:1:length(node_m)
            for i = 1:1:nB
                if any(strcmp(buildLocation{i,1},node_m))
                    build(end+1) = i;
                    if strcmp(networkNames{net},'Electrical')
                        Building(i).QPform.Electrical.subnetNode = n;
                        Building(i).QPform.Location = subNet.Electrical.Location(n);
                    elseif strcmp(networkNames{net},'DistrictHeat')
                        Building(i).QPform.DistrictHeat.subnetNode = n;
                        if ~isempty(subNet.DistrictHeat.Equipment{n})
                            equip = subNet.DistrictHeat.Equipment{n};
                            for j = 1:1:length(equip)
                                if ismember(Gen(equip(j)).Type,{'Heater';'CHP Generator';})
                                    Building(i).QPform.Heating = true;
                                end
                            end
                        end
                        if ~isempty(subNet.DistrictHeat.connections{n}) %connected to heaters at a different node
                            Building(i).QPform.Heating = true;
                        end
                        if Building(i).QPform.Heating
                            Building(i).QPform.H2E = Building(i).VariableStruct.FanPower/(1.025*(50-15))/1.225; %Flow = Heating/(Cp_Air*(Tsupply(t)-Tmix)); Power = FanPower*Flow/air density
                        end
                    elseif strcmp(networkNames{net},'DistrictCool')
                        Building(i).QPform.DistrictCool.subnetNode = n;    
                        if ~isempty(subNet.DistrictCool.Equipment{n})
                            equip = subNet.DistrictCool.Equipment{n};
                            for j = 1:1:length(equip)
                                if ismember(Gen(equip(j)).Type,{'Chiller';})
                                    Building(i).QPform.Cooling = true;
                                end
                            end
                        end
                        if ~isempty(subNet.DistrictHeat.connections{n}) %connected to heaters at a different node
                            Building(i).QPform.Cooling = true;
                        end
                        if Building(i).QPform.Cooling
                            Building(i).QPform.C2E = Building(i).VariableStruct.FanPower/(1.025*(25-12))/1.225; %Flow = Cooling/(Cp_Air*(Tmix-Tsupply(t))); Power = FanPower*Flow/air density
                        end
                    end
                end
            end
        end
        subNet.(networkNames{net}).Buildings{n} = build;
    end
end
end%Ends function locateBuildings