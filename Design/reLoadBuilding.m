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
    [InternalGains,~,Equipment,InteriorLighting, ExteriorLighting, OtherLoads] = BuildingLoads(Building(i),TestData.Weather.irradDireNorm,TestData.Weather.irradDiffHorz,Location,TestData.Timestamp);
    TestData.Building.NonHVACelectric(:,i) = Equipment + InteriorLighting + ExteriorLighting + OtherLoads;
    TestData.Building.InternalGains(:,i) = InternalGains;
end
