function loadBuilding
%% this loads a representation of a building a single zone temperature 
% Within the optimization 5 states are used:
% #1 represents the air zone temperature setpoint
% #2 represents the heating required to achieve that zone temperature setpoint
% #3 represents the cooling required to achieve that zone temperature setpoint
% #4 represents the temperature in excess of a comfort range (this is penalized)
% #5 represents the temperature below a comfort range (this is penalized)
global Plant
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
    for i = 1:1:nB
        f = .1;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
        UA_window = Plant.Building(i).VariableStruct.WallArea*Plant.Building(i).VariableStruct.WindowWallRatio*Plant.Building(i).VariableStruct.WindowUvalue/1000;
        UA_wall = (Plant.Building(i).VariableStruct.WallArea*(1-Plant.Building(i).VariableStruct.WindowWallRatio)/Plant.Building(i).VariableStruct.WallRvalue + Plant.Building(i).VariableStruct.RoofArea/Plant.Building(i).VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
        AirCapacitance = 2*Plant.Building(i).Area;
        
        QPform.states = {'T';'H';'C';'U';'L'};%Temperature, Heating, Cooling, exceeding upper comfort zone, exceeding lower comfort zone
        QPform.UA = UA_window+1/f*UA_wall;
        QPform.Cap = AirCapacitance;
        QPform.H2E = 1/Plant.Building(i).VariableStruct.COP_H;
        QPform.C2E = 1/Plant.Building(i).VariableStruct.COP_C;
        QPform.Cooling = false; %initial condition is that heating and cooling are converted to electricity, this can change in buildSubNet if there are local chillers or heaters
        QPform.Heating = false;
        QPform.Discomfort = Plant.Building(i).Area*Plant.Building(i).VariableStruct.occupancy; %cost of exceeding comfort band
        Plant.Building(i).QPform = QPform;
    end 
end