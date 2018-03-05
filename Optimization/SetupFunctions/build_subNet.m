function [subNet,Generator] = build_subNet(Network,Generator)
%identify generators and lines and their position in the network
%organize into sub-networks
% Group nodes between which transmission losses don't occur
% only create lines where transmission losses do occur
%% Network names, their corresponding abbreviation, what they represent
% Electrical --- E  --- standard 480V AC electrical network 
% DistrictHeat --- H --- standard 80C supply heating bus
% DistrictCool --- C --- standard 4C supply cooling bus
% Hydro        --- W --- River network with reservoirs and dams
% DirectCurrent --- DC --- 48V DC electrical network
% CoolingWater --- CW --- Water circulated between chillers and cooling towers
% Transmission1 --- E1 --- 230kV electric transmission (E2, E3, etc can be additional voltage levels
% Hydrogen     --- Hy --- Gaseous hydrogen stream
% LiqHydrogen  --- LH2 --- Liquid hydrogen
% Heating2     --- H2 --- Heat a different temperature than DistrictHeat (H3, H4... as needed)
%%-----%%%
nG = length(Generator);
nodes = length(Network);
genNames = cell(nG,1);
nodeNames = cell(nodes,1);
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nnAbrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
networkNames = fieldnames(Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));
for i = 1:1:nG
    genNames(i,1) = {Generator(i).Name};
end
listDams = {}; %%need better way to correlate columns of SourceSink data to dam
for i = 1:1:nodes
    nodeNames(i) = {Network(i).name};
    if isfield(Network,'Location') && ~isempty(Network(i).Location)
        location(i).Longitude = Network(i).Location.Longitude;
        location(i).Latitude = Network(i).Location.Latitude;
        location(i).TimeZone = Network(i).Location.TimeZone;
    else
        location(i) = {[]};
    end
    if isfield(Network,'Hydro') && ~isempty(Network(i).Hydro)
        listDams(end+1,1) = {Network(i).name};
    end
end
nLcum = 0; %cumulative line segment #
for net = 1:1:length(networkNames)
    subNet.(networkNames{net}).nodes = {};
    nIndex = nonzeros((1:length(nnList))'.*strcmp(networkNames{net},nnList));
    subNet.(networkNames{net}).abbreviation = nnAbrev{nIndex};
    subNet.(networkNames{net}).lineNames = {};
    subNet.(networkNames{net}).lineNumber = [];
    subNet.(networkNames{net}).lineLimit = [];
    if strcmp(networkNames{net},'Hydro')
        subNet.(networkNames{net}).lineMinimum = []; 
        subNet.(networkNames{net}).lineTime = []; 
    else
        subNet.(networkNames{net}).lineEff = [];
    end
    n = 0;
    for i = 1:1:nodes
        if ~isempty(Network(i).(networkNames{net}))
            if strcmp(networkNames{net},'Hydro')
                n = n+1;%add a new subnet node
                nLcum = nLcum+1;                
                subNet.Hydro.nodes(n) = {nodeNames(i)};
                subNet.Hydro.Location(n) = location(i);
                subNet.Hydro.nodeNumber(n) = nonzeros((1:length(listDams))'.*strcmp(nodeNames(i),listDams));%column index of this node in the stored matrices of Data.Hydro.SourceSink and Data.Hydro.Inflow
                subNet.Hydro.connections(n) = {Network(i).Hydro.connections};
                subNet.Hydro.Load(n) = {[]};
                subNet.Hydro.lineNumber(n,1) = nLcum;
                subNet.Hydro.lineMinimum(n,1) = Network(i).Hydro.InstreamFlow;
                subNet.Hydro.lineLimit(n,1) = inf;
                if ~isempty(Network(i).Hydro.connections)
                    subNet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_',Network(i).Hydro.connections);
                    subNet.Hydro.lineTime(n,1) = Network(i).Hydro.Time2Sea - Network(strcmp(Network(i).Hydro.connections,nodeNames)).Hydro.Time2Sea; %transit time from current river to downstream river
                else
                    subNet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_');%last dam before the sea
                    subNet.Hydro.lineTime(n,1) = Network(i).Hydro.Time2Sea; %no downstream dam
                end
            else
                %first check and see if this node is already part of a subNet node
                %nodes with perfect transmission are agregated into the first node in the nameList that they have perfect 2-way connection with
                [I,aNodes,connect] = agregatedNode(Network,nodeNames{i},networkNames{net});
                if I == i
                    n = n+1;%add a new subnet node
                    subNet.(networkNames{net}).nodes(n) = {aNodes};
                    subNet.(networkNames{net}).Location(n) = location(i);
                    L = [];
                    for j = 1:1:length(aNodes)
                        I = find(strcmp(aNodes{j},nodeNames),1,'first');
                        if isfield(Network(I).(networkNames{net}),'Load') && ~isempty(Network(I).(networkNames{net}).Load)%%note if there is a demand at this node
                            L(1,end+1) = Network(I).(networkNames{net}).Load;
                        end
                    end
                    subNet.(networkNames{net}).Load(n) = {L};
                    Connect = {};
                    for j=1:1:length(connect(:,1))
                        if ~any(strcmp(connect{j,2},aNodes))%imperfect transmission, need a line
                            [J, cNodes,~] = agregatedNode(Network,connect{j,2},networkNames{net});
                            pconnected = nodeNames{J};%name of node that the connected node will be agregated into if it is perfectly connected to any others
                            Connect(end+1) = {pconnected};
                            if J>i %new line connection, otherwise this was handled previously in the reverse direction
                                [eff, limit,dir] = lineProp(Network,subNet.(networkNames{net}).nodes{n},cNodes,networkNames{net});%find forward & reverse transmission efficiency & limit
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
    for n = 1:1:length(subNet.(networkNames{net}).nodes)
        gen = [];
        node_m = subNet.(networkNames{net}).nodes{n};
        for k = 1:1:length(node_m)
            i = nonzeros((1:length(nodeNames))'.*(strcmp(node_m{k},nodeNames)));
            equip = Network(i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                I = find(strcmp(equip{j}(s+1:end),genNames),1,'first');
                if ~isempty(I)
                    if isfield(Generator(I).QPform.output,subNet.(networkNames{net}).abbreviation)
                        gen(end+1) = I;
                        Generator(I).QPform.(networkNames{net}).subnetNode = n;
                    end
                else
                    disp(strcat('error, generator is not in library',equip{j}))
                end
            end
        end
        subNet.(networkNames{net}).Equipment{n} = gen;
    end
end
end%Ends function build SubNet

function [TransEff,TransLimit,dir] = lineProp(Network,node1,node2,net)
%find the transmission efficiency and limit between 2 connected nodes
%if one of the nodes is perfectly connected to another node there may be
%more than one pathway connecting them, so agregate the lines
nodes = length(Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Network(i).name};
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
        c = find(strcmp(node2{k},Network(I).(net).connections));
        if ~isempty(c)
            if isinf(Network(I).(net).Trans_Limit(c)) || TransLimit(1,1)==0
                TransEff(1,1) = Network(I).(net).Trans_Eff(c);
                TransLimit(1,1) = Network(I).(net).Trans_Limit(c);
            else
                weight = Network(I).(net).Trans_Limit(c)/(TransLimit(1,1) + Network(I).(net).Trans_Limit(c));
                TransEff(1,1) = (1-weight)*TransEff(1,1) + weight*Network(I).(net).Trans_Eff(c);%weighted efficiency of 2 concurrent lines
                TransLimit(1,1) = TransLimit(1,1) + Network(I).(net).Trans_Limit(c);
            end
        end
        %reverse direction efficiency and limit
        c = find(strcmp(node1{j},Network(J).(net).connections));
        if ~isempty(c)
            if isinf(Network(J).(net).Trans_Limit(c)) || TransLimit(1,2)==0
                TransEff(1,2) = Network(J).(net).Trans_Eff(c);
                TransLimit(1,2) = Network(J).(net).Trans_Limit(c);
            else
                weight = Network(J).(net).Trans_Limit(c)/(TransLimit(1,2) + Network(J).(net).Trans_Limit(c));
                TransEff(1,2) = (1-weight)*TransEff(1,2) + weight*Network(J).(net).Trans_Eff(c);
                TransLimit(1,2) = TransLimit(1,2) + Network(J).(net).Trans_Limit(c);
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

function [I,aNodes,connect] = agregatedNode(Network,node,net)
%Any connected nodes with perfect bi-directional transfer are agregated into the node earliest in the list of names
%This function finds which node in the list that is
nodes = length(Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Network(i).name};
end
I = find(strcmp(node,nodeNames),1,'first');
aNodes = {node};
connect = cell(length(Network(I).(net).connections),2);
connect(:,2) = Network(I).(net).connections;
connect(:,1) = {node};
[m,~] = size(connect);
k = 0;
while k<m
    k = k+1;
    if ~any(strcmp(connect(k,2),aNodes))%avoid looking at connections to nodes already in agregated node
        [TransEff, ~,~] = lineProp(Network,connect(k,1),connect(k,2),net);
        if ~isempty(TransEff) && length(TransEff) == 2 && min(TransEff)==1 && ~strcmp(net,'Hydro')%perfect bi-directional energy transfer, hydro lines are river segments, can't agregate
            J = find(strcmp(connect{k,2},nodeNames),1,'first');
            aNodes(end+1) = nodeNames(J);%add to list of agregated nodes
            %%add additional connections to check
            c = length(Network(J).(net).connections);
            connect(end+1:end+c,1) = connect(k,2);
            connect(end-c+1:end,2) = Network(J).(net).connections;
            I = min(J,I); %keep lowest number (index in list of node names
        end
        [m,~] = size(connect);
    end
end
end%Ends function agregateNode