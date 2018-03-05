function reLoadBuilding(Building,Network)
global TestData
nB = length(Building);
nS = length(TestData.Timestamp);
TestData.Building.NonHVACelectric = zeros(nS,nB);
TestData.Building.InternalGains = zeros(nS,nB);
nodes = length(Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Network(i).name};
end
for i = 1:1:nB
    I = nonzeros((1:nodes)'.*(strcmp(Building(i).Location,nodeNames)));
    if isfield(Network(I),'Location')
        Location = Network(I).Location;
    else
        Location = {[]};
    end
    SG = SolarGain(Building(i),TestData.Timestamp,Location,TestData.Weather);
    B_loads = BuildingLoads(Building(i),TestData.Timestamp,SG);
    TestData.Building.InternalGains(:,i) = B_loads.InternalGains;
    TestData.Building.NonHVACelectric(:,i) = B_loads.Equipment + B_loads.InteriorLighting + B_loads.ExteriorLighting + B_loads.OtherLoads;
    if ~isempty(B_loads.DCloads)
        TestData.Building.DCloads(:,i) = B_loads.DCloads;
    end
end
