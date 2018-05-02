function [subnet,gen] = load_network(network,gen)
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
nG = length(gen);
nodes = length(network);
genNames = cell(nG,1);
nodeNames = cell(nodes,1);
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nnAbrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
networkNames = fieldnames(network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));
for i = 1:1:nG
    genNames(i,1) = {gen(i).Name};
end
listDams = {}; %%need better way to correlate columns of SourceSink data to dam
for i = 1:1:nodes
    nodeNames(i) = {network(i).name};
    if isfield(network,'Location') && ~isempty(network(i).Location)
        location(i).Longitude = network(i).Location.Longitude;
        location(i).Latitude = network(i).Location.Latitude;
        location(i).TimeZone = network(i).Location.TimeZone;
    else
        location(i) = {[]};
    end
    if isfield(network,'Hydro') && ~isempty(network(i).Hydro)
        listDams(end+1,1) = {network(i).name};
    end
end
line_number = 0; %cumulative line segment #
for net = 1:1:length(networkNames)
    subnet.(networkNames{net}).nodes = {};
    n_index = nonzeros((1:length(nnList))'.*strcmp(networkNames{net},nnList));
    subnet.(networkNames{net}).abbreviation = nnAbrev{n_index};
    subnet.(networkNames{net}).lineNames = {};
    subnet.(networkNames{net}).lineNumber = [];
    subnet.(networkNames{net}).lineLimit = [];
    if strcmp(networkNames{net},'Hydro')
        subnet.(networkNames{net}).lineMinimum = []; 
        subnet.(networkNames{net}).lineTime = []; 
    else
        subnet.(networkNames{net}).lineEff = [];
    end
    n = 0;
    for i = 1:1:nodes
        if ~isempty(network(i).(networkNames{net}))
            if strcmp(networkNames{net},'Hydro')
                n = n+1;%add a new subnet node
                line_number = line_number+1;                
                subnet.Hydro.nodes(n) = {nodeNames(i)};
                subnet.Hydro.Location(n) = location(i);
                subnet.Hydro.nodeNumber(n) = nonzeros((1:length(listDams))'.*strcmp(nodeNames(i),listDams));%column index of this node in the stored matrices of Data.Hydro.SourceSink and Data.Hydro.Inflow
                subnet.Hydro.connections(n) = {network(i).Hydro.connections};
                subnet.Hydro.Load(n) = {[]};
                subnet.Hydro.lineNumber(n,1) = line_number;
                subnet.Hydro.lineMinimum(n,1) = network(i).Hydro.InstreamFlow;
                subnet.Hydro.lineLimit(n,1) = inf;
                if ~isempty(network(i).Hydro.connections)
                    subnet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_',network(i).Hydro.connections);
                    subnet.Hydro.lineTime(n,1) = network(i).Hydro.Time2Sea - network(strcmp(network(i).Hydro.connections,nodeNames)).Hydro.Time2Sea; %transit time from current river to downstream river
                else
                    subnet.Hydro.lineNames(n,1) = strcat(nodeNames(i),'_Hydro_');%last dam before the sea
                    subnet.Hydro.lineTime(n,1) = network(i).Hydro.Time2Sea; %no downstream dam
                end
            else
                %first check and see if this node is already part of a subNet node
                %nodes with perfect transmission are agregated into the first node in the nameList that they have perfect 2-way connection with
                [I,aNodes,connect] = agregated_node(network,nodeNames{i},networkNames{net});
                if I == i
                    n = n+1;%add a new subnet node
                    subnet.(networkNames{net}).nodes(n) = {aNodes};
                    subnet.(networkNames{net}).Location(n) = location(i);
                    L = [];
                    for j = 1:1:length(aNodes)
                        I = find(strcmp(aNodes{j},nodeNames),1,'first');
                        if isfield(network(I).(networkNames{net}),'Load') && ~isempty(network(I).(networkNames{net}).Load)%%note if there is a demand at this node
                            L(1,end+1) = network(I).(networkNames{net}).Load;
                        end
                    end
                    subnet.(networkNames{net}).Load(n) = {L};
                    Connect = {};
                    for j=1:1:length(connect(:,1))
                        if ~any(strcmp(connect{j,2},aNodes))%imperfect transmission, need a line
                            [J, cNodes,~] = agregated_node(network,connect{j,2},networkNames{net});
                            pconnected = nodeNames{J};%name of node that the connected node will be agregated into if it is perfectly connected to any others
                            Connect(end+1) = {pconnected};
                            if J>i %new line connection, otherwise this was handled previously in the reverse direction
                                [eff, limit,dir] = line_prop(network,subnet.(networkNames{net}).nodes{n},cNodes,networkNames{net});%find forward & reverse transmission efficiency & limit
                                if strcmp(dir,'none') %no transmission (zero efficiency)
                                    %do nothing
                                else
                                    line_number = line_number+1;
                                    if strcmp(dir,'reverse')
                                        subnet.(networkNames{net}).lineNames(end+1,1) = (strcat(pconnected,'_',networkNames{net},'_',nodeNames(i)));
                                    else
                                        subnet.(networkNames{net}).lineNames(end+1,1) = (strcat(nodeNames(i),'_',networkNames{net},'_',pconnected));
                                    end
                                    if strcmp(dir,'dual')
                                        subnet.(networkNames{net}).lineEff(end+1,1:2) = eff;
                                        subnet.(networkNames{net}).lineLimit(end+1,1:2) = limit;
                                    else
                                        subnet.(networkNames{net}).lineEff(end+1,1) = eff;
                                        subnet.(networkNames{net}).lineLimit(end+1,1) = limit;
                                    end
                                    subnet.(networkNames{net}).lineNumber(end+1,1) = line_number;
                                end
                            end
                        end
                    end
                    subnet.(networkNames{net}).connections{n} = Connect;
                end
            end
        end
    end
end

for net = 1:1:length(networkNames)
%identify equipment at each subNet node (equipment can apear in multiple
%sub-nets if it produces heat and power, or uses water to produce electricity
    for n = 1:1:length(subnet.(networkNames{net}).nodes)
        gen_at_node = [];
        node_m = subnet.(networkNames{net}).nodes{n};
        for k = 1:1:length(node_m)
            i = nonzeros((1:length(nodeNames))'.*(strcmp(node_m{k},nodeNames)));
            equip = network(i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                I = find(strcmp(equip{j}(s+1:end),genNames),1,'first');
                if ~isempty(I)
                    if isfield(gen(I).QPform.output,subnet.(networkNames{net}).abbreviation)
                        gen_at_node(end+1) = I;
                        gen(I).QPform.(networkNames{net}).subnetNode = n;
                    end
                else
                    disp(strcat('error, generator is not in library',equip{j}))
                end
            end
        end
        subnet.(networkNames{net}).Equipment{n} = gen_at_node;
    end
end
end%Ends function load_network

function [TransEff,TransLimit,dir] = line_prop(Network,node1,node2,net)
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

function [I,aNodes,connect] = agregated_node(Network,node,net)
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
        [TransEff, ~,~] = line_prop(Network,connect(k,1),connect(k,2),net);
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