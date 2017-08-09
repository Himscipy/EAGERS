function Power = findHydroPower(gen)
%given the generator index, find the downstream flow and spill flow to
%calculate the power
global Plant DateSim
upriver = {};
spill = {};
upLines = [];
spillLines = [];
for i = 1:1:length(Plant.subNet.lineNames.Hydro)
    name = Plant.subNet.lineNames.Hydro{i};
    k = strfind(name,'_');
    if strcmp(name(k(1)+1:k(2)-1),'Spill')
        spill(end+1) = {name(1:k(1)-1)};
        spillLines(end+1) = i;
    else
        upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
        upLines(end+1) = i;
    end
end
for i = 1:1:length(Plant.subNet.Hydro)
    equip = Plant.subNet.Hydro(i).Equipment;
    if any(equip==gen) %this is the node wih this generator
        nodeName = Plant.subNet.Hydro(i).nodes{1};
    end
end
J = upLines(strcmp(nodeName,upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
S = spillLines(strcmp(nodeName,spill));% spill flow of this node, subtracts from the power generation  (should be at most 1)
OutFlow = getHydroFlows(DateSim,J);
SpillFlow = getHydroFlows(DateSim,S);
Power = (OutFlow - SpillFlow)*Plant.Generator(gen).OpMatA.output.E;%Converting water flow to power            