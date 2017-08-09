function build_subNet
%identify generators and lines and their position in the network
%organize into sub-networks
% Group nodes between which transmission losses don't occur
% only create lines where transmission losses do occur
global Plant 
nG = length(Plant.Generator);
nB = length(Plant.Building);
nodes = length(Plant.Network);
genNames = cell(nG,1);
nodeNames = cell(nodes,1);
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));

for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
    nodeDirectory(i).Name = Plant.Network(i).name;
    if isfield(Plant.Network,'Location') && ~isempty(Plant.Network(i).Location)
        location(i).Longitude = Plant.Network(i).Location.Longitude;
        location(i).Latitude = Plant.Network(i).Location.Latitude;
        location(i).TimeZone = Plant.Network(i).Location.TimeZone;
    else location(i) = {[]};
    end
end
nLcum = 0; %cumulative line segment #
for net = 1:1:length(networkNames)
    subNet.(networkNames{net}).nodes = {};
    subNet.(networkNames{net}).lineNames = {};
    subNet.(networkNames{net}).lineNumber = [];
    subNet.(networkNames{net}).lineLimit = [];
    if strcmp(networkNames{net},'Hydro')
        subNet.(networkNames{net}).lineMinimum = []; 
        subNet.(networkNames{net}).lineTime = []; 
    else
        subNet.(networkNames{net}).lineEff = [];   %added space for line effeciencies
    end
    n = 0;
    for i = 1:1:nodes
        if ~isempty(Plant.Network(i).(networkNames{net}))
            if strcmp(networkNames{net},'Hydro')
                n = n+1;%add a new subnet node
                nLcum = nLcum+1;                
                subNet.Hydro.nodes(n) = {nodeNames(i)};
                subNet.Hydro.Location(n) = location(i);
                subNet.Hydro.nodeNumber(n) = find(strcmp(nodeNames(i),Plant.Data.Hydro.Nodes));%column index of this node in the stored matrices of Data.Hydro.SourceSink and Data.Hydro.Inflow
                subNet.Hydro.connections(n) = {Plant.Network(i).Hydro.connections};
                subNet.Hydro.Load(n) = {[]};
                subNet.Hydro.lineNumber(n,1) = nLcum;
                subNet.Hydro.lineMinimum(n,1) = Plant.Network(i).Hydro.InstreamFlow;
                subNet.Hydro.lineLimit(n,1) = inf;
                if ~isempty(Plant.Network(i).Hydro.connections)
                    subNet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_',Plant.Network(i).Hydro.connections);
                    subNet.Hydro.lineTime(n,1) = Plant.Network(i).Hydro.Time2Sea - Plant.Network(strcmp(Plant.Network(i).Hydro.connections,nodeNames)).Hydro.Time2Sea; %transit time from current river to downstream river
                else
                    subNet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_');%last dam before the sea
                    subNet.Hydro.lineTime(n,1) = Plant.Network(i).Hydro.Time2Sea; %no downstream dam
                end
            else
                %first check and see if this node is already part of a subNet node
                %nodes with perfect transmission are agregated into the first node in the nameList that they have perfect 2-way connection with
                [I,aNodes,connect] = agregatedNode(nodeNames{i},networkNames{net});
                if I == i
                    n = n+1;%add a new subnet node
                    subNet.(networkNames{net}).nodes(n) = {aNodes};
                    subNet.(networkNames{net}).Location(n) = location(i);
                    L = [];
                    for j = 1:1:length(aNodes)
                        I = find(strcmp(aNodes{j},nodeNames),1,'first');
%                         L =[]; %This deletes the stored load value before storing it; Moved to line 71
                        if isfield(Plant.Network(I).(networkNames{net}),'Load') && ~isempty(Plant.Network(I).(networkNames{net}).Load)%%note if there is a demand at this node
                            L(1,end+1) = Plant.Network(I).(networkNames{net}).Load;
                        end
                    end
                    subNet.(networkNames{net}).Load(n) = {L};
                    Connect = {};
                    for j=1:1:length(connect(:,1))
                        if ~any(strcmp(connect{j,2},aNodes))%imperfect transmission, need a line
                            [J, cNodes,~] = agregatedNode(connect{j,2},networkNames{net});
                            pconnected = nodeNames{J};%name of node that the connected node will be agregated into if it is perfectly connected to any others
                            Connect(end+1) = {pconnected};
                            if J>i %new line connection, otherwise this was handled previously in the reverse direction
                                [eff, limit,dir] = lineProp(subNet.(networkNames{net}).nodes{n},cNodes,networkNames{net});%find forward & reverse transmission efficiency & limit
                                if strcmp(dir,'none') %no transmission (zero efficiency)
                                    %do nothing
                                else
                                    nLcum = nLcum+1;
                                    if strcmp(dir,'reverse')
                                        subNet.(networkNames{net}).lineNames(end+1,1) = (strcat(pconnected,'_',networkNames{net},'_',nodeNames(i)));
                                    else
                                        subNet.(networkNames{net}).lineNames(end+1,1) = (strcat(nodeNames(i),'_',networkNames{net},'_',pconnected));
                                    end
                                    if strcmp(dir,'dual')
                                        subNet.(networkNames{net}).lineEff(end+1,1:2) = eff;
                                        subNet.(networkNames{net}).lineLimit(end+1,1:2) = limit;
                                    else
                                        subNet.(networkNames{net}).lineEff(end+1,1) = eff;
                                        subNet.(networkNames{net}).lineLimit(end+1,1) = limit;
                                    end
                                    subNet.(networkNames{net}).lineNumber(end+1,1) = nLcum;
                                end
                            end
                        end
                    end
                    subNet.(networkNames{net}).connections{n} = Connect;
                end
            end
        end
    end
end

for net = 1:1:length(networkNames)
%identify equipment at each subNet node (equipment can apear in multiple
%sub-nets if it produces heat and power, or uses water to produce electricity
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
    elseif strcmp(networkNames{net},'Hydro')
        out = 'W';
    end
    for n = 1:1:length(subNet.(networkNames{net}).nodes)
        gen = [];
        node_m = subNet.(networkNames{net}).nodes{n};
        for k = 1:1:length(node_m)
            i = find(strcmp(node_m{k},nodeNames),1,'first');
            nodeDirectory(i).(networkNames{net}) = n;
            equip = Plant.Network(i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                I = find(strcmp(equip{j}(s+1:end),genNames),1,'first');
                if ~isempty(I)
                    if isfield(Plant.Generator(I).QPform.output,out)
                        gen(end+1) = I;
                        Plant.Generator(I).QPform.(networkNames{net}).subnetNode = n;
                    end
                else disp(strcat('error, generator is not in library',equip{j}))
                end
            end
        end
        subNet.(networkNames{net}).Equipment{n} = gen;
    end
end
Plant.subNet = subNet;

%identify upper bound for utility states
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Utility') && ~isempty(Plant.Generator(i).QPform.states)%%avoid things like gas utility with no states
        %identify the network
        out = char(fieldnames(Plant.Generator(i).QPform.output));
        if strcmp(out,'E')
            net = 'Electrical';
        elseif strcmp(out,'H')
            net = 'DistrictHeat';
        elseif strcmp(out,'C')
            net = 'DistrictCool';
        elseif strcmp(out,'W')
            net = 'Hydro';
        end
        if isfield(Plant.Data,'Demand')
            Plant.Generator(i).QPform.X.ub = max(Plant.Data.Demand.(out));%max Purchase
        else
            Plant.Generator(i).QPform.X.ub = 1e10; %arbitrary upper bound that is not inf
        end
        if isfield(Plant.Generator(i).QPform,'Y')
            maxSellback = 0;
            for n = 1:1:length(subNet.(net).nodes)
                equip = subNet.(net).Equipment{n};
                for j = 1:1:length(equip)
                    maxSellback = maxSellback + max(0,Plant.Generator(equip(j)).Size*Plant.Generator(equip(j)).QPform.output.(out)(1));
                end
            end
            Plant.Generator(i).QPform.Y.ub = maxSellback;
        end
    end
end  

%identify the location of any buildings & if they are connected to heaters and chillers
for i = 1:1:nB
    I = nonzeros((1:nB)'.*(strcmp(Plant.Building(i).Location,nodeNames)));
    Plant.Building(i).QPform.nodeE = nodeDirectory(I).Electrical; %node in the electric network
    Plant.Building(i).QPform.Location = subNet.Electrical.Location(Plant.Building(i).QPform.nodeE);
    if isfield(subNet,'DistrictHeat')
        Plant.Building(i).QPform.nodeH = nodeDirectory(I).DistrictHeat; %node in the electric network
        if ~isempty(subNet.DistrictHeat.connections{Plant.Building(i).QPform.nodeH}) %connected to heaters at a different node
            Plant.Building(i).QPform.Heating = true;
        else
            equip = subNet.DistrictHeat.Equipment{Plant.Building(i).QPform.nodeH};
            for k = 1:1:length(equip)
                if ismember(Plant.Generator(equip(k)).Type,{'Heater';'CHP Generator';})
                    Plant.Building(i).QPform.Heating = true;
                end
            end
        end
        if Plant.Building(i).QPform.Heating
            %%Need to update Plant.Building(i).QPform.H2E to only be
            %%the fans and other associated electric loads, rather than
            %%the COP of heating
        end
    end
    if isfield(subNet,'DistrictCool')
        Plant.Building(i).QPform.nodeC = nodeDirectory(I).DistrictCool; %node in the electric network
        if ~isempty(subNet.DistrictHeat.connections{Plant.Building(i).QPform.nodeC}) %connected to heaters at a different node
            Plant.Building(i).QPform.Cooling = true;
        else
            equip = subNet.DistrictCool.Equipment{Plant.Building(i).QPform.nodeC};
            for k = 1:1:length(equip)
                if ismember(Plant.Generator(equip(k)).Type,{'Chiller';})
                    Plant.Building(i).QPform.Cooling = true;
                end
            end
        end
        if Plant.Building(i).QPform.Cooling
            %%Need to update Plant.Building(i).QPform.C2E to only be
            %%the fans and other associated electric loads, rather than
            %%the COP of cooling
        end
    end    
end

end%Ends function build SubNet

function [TransEff,TransLimit,dir] = lineProp(node1,node2,net)
%find the transmission efficiency and limit between 2 connected nodes
%if one of the nodes is perfectly connected to another node there may be
%more than one pathway connecting them, so agregate the lines
global Plant
nodes = length(Plant.Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
TransEff = zeros(1,2);
TransLimit = zeros(1,2);
node1 = unique(node1);
node2 = unique(node2);
for j = 1:1:length(node1)
    I = find(strcmp(node1{j},nodeNames),1,'first');
    for k = 1:1:length(node2)
        J = find(strcmp(node2{k},nodeNames),1,'first');
        %forward direction efficieny & limit
        c = find(strcmp(node2{k},Plant.Network(I).(net).connections));
        if ~isempty(c)
            if isinf(Plant.Network(I).(net).Trans_Limit(c)) || TransLimit(1,1)==0
                TransEff(1,1) = Plant.Network(I).(net).Trans_Eff(c);
                TransLimit(1,1) = Plant.Network(I).(net).Trans_Limit(c);
            else
                weight = Plant.Network(I).(net).Trans_Limit(c)/(TransLimit(1,1) + Plant.Network(I).(net).Trans_Limit(c));
                TransEff(1,1) = (1-weight)*TransEff(1,1) + weight*Plant.Network(I).(net).Trans_Eff(c);%weighted efficiency of 2 concurrent lines
                TransLimit(1,1) = TransLimit(1,1) + Plant.Network(I).(net).Trans_Limit(c);
            end
        end
        %reverse direction efficiency and limit
        c = find(strcmp(node1{j},Plant.Network(J).(net).connections));
        if ~isempty(c)
            if isinf(Plant.Network(J).(net).Trans_Limit(c)) || TransLimit(1,2)==0
                TransEff(1,2) = Plant.Network(J).(net).Trans_Eff(c);
                TransLimit(1,2) = Plant.Network(J).(net).Trans_Limit(c);
            else
                weight = Plant.Network(J).(net).Trans_Limit(c)/(TransLimit(1,2) + Plant.Network(J).(net).Trans_Limit(c));
                TransEff(1,2) = (1-weight)*TransEff(1,2) + weight*Plant.Network(J).(net).Trans_Eff(c);
                TransLimit(1,2) = TransLimit(1,2) + Plant.Network(J).(net).Trans_Limit(c);
            end
        end
    end
end 
if TransEff(1,1)==0 && TransEff(1,2)>0
    dir = 'reverse';
    TransEff = TransEff(1,2);
    TransLimit = TransLimit(1,2);
elseif TransEff(1,1)>0 && TransEff(1,2)==0
    dir = 'forward';
    TransEff = TransEff(1,1);
    TransLimit = TransLimit(1,1);
elseif TransEff(1,1)==0 && TransEff(1,2)==0
    dir = 'none';
    TransEff = [];
    TransLimit = [];
else
    dir = 'dual';
end
end%Ends function lineProp

function [I,aNodes,connect] = agregatedNode(node,net)
%Any connected nodes with perfect bi-directional transfer are agregated into the node earliest in the list of names
%This function finds which node in the list that is
global Plant
nodes = length(Plant.Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
I = find(strcmp(node,nodeNames),1,'first');
aNodes = {node};
connect = cell(length(Plant.Network(I).(net).connections),2);
connect(:,2) = Plant.Network(I).(net).connections;
connect(:,1) = {node};
[m,~] = size(connect);
k = 0;
while k<m
    k = k+1;
    if ~any(strcmp(connect(k,2),aNodes))%avoid looking at connections to nodes already in agregated node
        [TransEff, ~,~] = lineProp(connect(k,1),connect(k,2),net);
        if ~isempty(TransEff) && length(TransEff) == 2 && min(TransEff)==1 && ~strcmp(net,'Hydro')%perfect bi-directional energy transfer, hydro lines are river segments, can't agregate
            J = find(strcmp(connect{k,2},nodeNames),1,'first');
            aNodes(end+1) = nodeNames(J);%add to list of agregated nodes
            %%add additional connections to check
            c = length(Plant.Network(J).(net).connections);
            connect(end+1:end+c,1) = connect(k,2);
            connect(end-c+1:end,2) = Plant.Network(J).(net).connections;
            I = min(J,I); %keep lowest number (index in list of node names
        end
        [m,~] = size(connect);
    end
end
end%Ends function agregateNode