function loadBuilding% Loads generators for economic dispatch
%% this loads a representation of a building with 6 states
global Plant
nB = length(Plant.Building);
for i = 1:1:nB
    QPform.states = {'T';'H';'C';'U';'L'};%Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
    QPform.UA = Plant.Building(i).Area/Plant.Building(i).VariableStruct.Resistance;
    QPform.Cap = Plant.Building(i).Area*Plant.Building(i).VariableStruct.Capacitance;
    QPform.H2E = 1/Plant.Building(i).VariableStruct.COP_H;
    QPform.C2E = 1/Plant.Building(i).VariableStruct.COP_C;
    QPform.Cooling = false; %initial condition is that heating and cooling are converted to electricity, this can change in buildSubNet if there are local chillers or heaters
    QPform.Heating = false;
    QPform.Discomfort = Plant.Building(i).Area*Plant.Building(i).VariableStruct.occupancy; %cost of exceeding comfort band
    Plant.Building(i).QPform = QPform;
end 