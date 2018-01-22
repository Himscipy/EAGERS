global Plant Model_dir
Plant.Name = 'ColumbiaRiverDams';
%%Columbia river: 'GC';'CJ';'W';'RR';'RI';'WNPM';'PR';'MN';'JD';'TD';'BNVL';
%%Snake River: 'BRWNL';'OX';'HC';'LWRGRNT';'LTLGSE';'LWRMONUM';'ICEHRBR';

%%rows 1->366 (hourly): October 2nd, 2015 -> October 1st, 2016; 
Plant.Data.Timestamp = linspace(datenum([2015 10 2 1 0 0]),datenum([2016 10 2 0 0 0]),8784);
Plant.Data.Weather.Tdb = 20*ones(length(Plant.Data.Timestamp),1);

%% Columns of SourceSink, SpillFlow, InFlow and Outflow must correspond to the node in the order of NodeNames
%%columns 1->18 (kcfs or KW): Grand Coulee, Chief Joseph, Wells, Rocky Ridge, Rock Island, Wanapum, PriestRiver, McNary, John Day, The Dalles, Bonneville, Brownlee, Oxbow, Hells Canyon, Lower Granite, Little Goose, Lower Monumental, Ice Harbor;
%SourceSink = Mass balance; Sinks and Sources in river segments; Negative values = sinks; Positive values = sources;
load(fullfile(Model_dir,'Data','Hydro','SourcesandSinks')); load(fullfile(Model_dir,'Data','Hydro','Powerflow'));load(fullfile(Model_dir,'Data','Hydro','Spillflow'));load(fullfile(Model_dir,'Data','Hydro','Outflow'));load(fullfile(Model_dir,'Data','Hydro','Inflow'));load(fullfile(Model_dir,'Data','Hydro','PowerGen'));
Plant.Data.Hydro.SourceSink = SourcesandSinks; 
%Powerflow = Flow used to produce power
Plant.Data.Hydro.PowerFlow = Powerflow;
%Spillflow = Flow that is spilled not used for power
Plant.Data.Hydro.SpillFlow = Spillflow;
%Outflow = Full amount of discharge flow
Plant.Data.Hydro.OutFlow = Outflow;
%columns (kilo-cfs): inflow(GC), infl(CJ)-Disch(GC), infl(W)-Disch(CJ), infl(RR)-Disch(W), infl(RI)-Disch(RR), infl(WNP)-Disch(RI), infl(PR)-Disch(WNP), infl(MN)-Disch(IH)-Disch(PR), infl(JD)-Disch(MN), infl(TD)-Disch(JD), infl(BNVL)-Disch(TD), inflow(BRW), infl(OX)-Disch(BRW), infl(HC)-Disch(OX), infl(LGRT)-Disch(HC), infl(LGS)-Disch(LGRT), infl(LMNT)-Disch(LG), infl(IH)-Disch(LMNT);
Plant.Data.Hydro.InFlow = Inflow;
%PowerGen = Generation (kW) every hour
Plant.Data.Hydro.PowerGen = PowerGen;
PowerGen(isnan(PowerGen)) = 0;
Plant.Data.Demand = [];
%Still need 
%            storage data

NodeNames = {'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa';'Wenatchee, Wa';...
             'South Wenatchee, Wa';'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';...
             'Rufus, Or';'The Dalles, Or';'Bonneville, Or'; 'Cambridge, Id';...
             'Oxbow, Or'; 'Hells Canyon, Or'; 'Almota, Wa'; 'Starbuck, Wa';...
             'Kahlotus, Wa';'Pasco, Wa';}; %18 nodes in the network (currently overlapping electrical and hydro)
         
DownRiver = {'Bridgeport, Wa';'Azwell, Wa';'Wenatchee, Wa';'South Wenatchee, Wa';...
             'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';'Rufus, Or';...
             'The Dalles, Or';'Bonneville, Or'; '';'Oxbow, Or'; 'Hells Canyon, Or';...
             'Almota, Wa'; 'Starbuck, Wa';'Kahlotus, Wa';'Pasco, Wa';'Umatilla, Or';}; %Water connections (& electrical connections)
         
UpRiver = {{''};{'Grand Coulee, Wa'};{'Bridgeport, Wa'};{'Azwell, Wa'};{'Wenatchee, Wa'};...
             {'South Wenatchee, Wa'};{'Beverly, Wa'};{'Mattawa, Wa','Pasco, Wa'};...
             {'Umatilla, Or'};{'Rufus, Or'};{'The Dalles, Or'};{''};{'Cambridge, Id'};...
             {'Oxbow, Or'}; {'Hells Canyon, Or'};{'Almota, Wa'};{'Starbuck, Wa'};...
             {'Kahlotus, Wa'};}; %Electrical network connections

Equipment = {'Grand Coulee';'Chief Joseph';'Wells';'Rocky Reach';'Rock Island';...
             'Wanapum';'Priest Rapids';'McNary';'John Day';'The Dalles';...
             'Bonneville';'Brownlee';'Oxbow';'Hells Canyon';'Lower Granite';...
             'Little Goose';'Lower Monumental';'Ice Harbor';}; %Equipment at each node

InstreamFlow = [ 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;]; %instream flow requirements at these 18 node points

Time2Sea = [80.6;73.7;69.7;64.0;61.3;56.2;53.7;39.5;29.1;25.9;19.7;82.4;80.8;77.4;58.5;53.4;49.5;45.2];%Hours for water to reach mouth of columbia from each node

RiverMile = [596.6;545.1;515.8;473.7;453.4;415.8;397.1;292;215.6;191.5;146.1;285;273;247.6;107.5;70.3;41.6;9.7;]; %was used to calculate the Time2Sea

%optimization options
Plant.optimoptions.Interval = 31;
Plant.optimoptions.Horizon = 24*7;
Plant.optimoptions.Resolution = 4;
Plant.optimoptions.Topt = 600;
Plant.optimoptions.Tmpc = 60;
Plant.optimoptions.nsSmooth = 0;
Plant.optimoptions.scaletime = 1;
Plant.optimoptions.tspacing = 'constant';
Plant.optimoptions.sequential = false;
Plant.optimoptions.excessHeat = true;
Plant.optimoptions.thresholdSteps = 1;
Plant.optimoptions.Buffer = 40;
Plant.optimoptions.method = 'Dispatch';
Plant.optimoptions.MixedInteger = false;
Plant.optimoptions.SpinReserve = false;
Plant.optimoptions.SpinReservePerc = 0;
Plant.optimoptions.forecast = 'Perfect';
Plant.optimoptions.solver = 'quadprog';
Plant.optimoptions.mode = 'virtual';

Dam(1).Type = 'Hydro Storage';
Dam(1).Name = 'Grand Coulee';
Dam(1).Source = 'Water';
Dam(1).Output = [];
Dam(1).Size = 9170;%This would be storage capacity in kilo-acre-feet
Dam(1).Enabled = 1;
Dam(1).VariableStruct.MaxGenCapacity = 4414000; %power in kW
Dam(1).VariableStruct.RampUp = 189; %power change in kW/hr
Dam(1).VariableStruct.RampDown = 233;
Dam(1).VariableStruct.MaxGenFlow = 202.5; %flow in 1000 cfs
Dam(1).VariableStruct.MaxSpillFlow = 0.2; %flow in 1000 cfs
Dam(1).VariableStruct.MaxHead = 334.2;
Dam(1).VariableStruct.MinHead = 180; %guess minimum height allowable

Dam(2).Type = 'Hydro Storage';
Dam(2).Name = 'Chief Joseph'; 
Dam(2).Source = 'Water';
Dam(2).Output = [];
Dam(2).Size = 588; %This would be storage capacity in kilo-acre-feet
Dam(2).Enabled = 1;
Dam(2).VariableStruct.MaxGenCapacity = 2180000; %power in kW
Dam(2).VariableStruct.RampUp =532; %power change in kW/hr
Dam(2).VariableStruct.RampDown =586; 
Dam(2).VariableStruct.MaxGenFlow = 179.3; %flow in 1000 cfs
Dam(2).VariableStruct.MaxSpillFlow = 31.3; %flow in 1000 cfs
Dam(2).VariableStruct.MaxHead = 180.8; 
Dam(2).VariableStruct.MinHead = 100; %guess minimum height allowable

Dam(3).Type = 'Hydro Storage';
Dam(3).Name = 'Wells';
Dam(3).Source = 'Water';
Dam(3).Output = [];
Dam(3).Size = 98; %This would be storage capacity in kilo-acre-feet
Dam(3).Enabled = 1;
Dam(3).VariableStruct.MaxGenCapacity = 742000; %power in kW
Dam(3).VariableStruct.RampUp =533; %power change in kW/hr
Dam(3).VariableStruct.RampDown =608; 
Dam(3).VariableStruct.MaxGenFlow = 175.3; %flow in 1000 cfs
Dam(3).VariableStruct.MaxSpillFlow = 57.96; %flow in 1000 cfs
Dam(3).VariableStruct.MaxHead = 75.19; 
Dam(3).VariableStruct.MinHead = 40; %guess minimum height allowable

Dam(4).Type = 'Hydro Storage';
Dam(4).Name =  'Rocky Reach'; 
Dam(4).Source = 'Water';
Dam(4).Output = [];
Dam(4).Size = 37; %This would be storage capacity in kilo-acre-feet
Dam(4).Enabled = 1;
Dam(4).VariableStruct.MaxGenCapacity = 1012000; %power in kW
Dam(4).VariableStruct.RampUp =396; %power change in kW/hr
Dam(4).VariableStruct.RampDown =476; 
Dam(4).VariableStruct.MaxGenFlow = 159.8; %flow in 1000 cfs
Dam(4).VariableStruct.MaxSpillFlow = 52.99; %flow in 1000 cfs
Dam(4).VariableStruct.MaxHead =  95.19; 
Dam(4).VariableStruct.MinHead = 60; %guess minimum height allowable
 
Dam(5).Type = 'Hydro Storage';
Dam(5).Name = 'Rock Island'; 
Dam(5).Source = 'Water';
Dam(5).Output = [];
Dam(5).Size = 12; %This would be storage capacity in kilo-acre-feet
Dam(5).Enabled = 1;
Dam(5).VariableStruct.MaxGenCapacity = 449000; %power in kW
Dam(5).VariableStruct.RampUp = 391; %power change in kW/hr
Dam(5).VariableStruct.RampDown = 478; 
Dam(5).VariableStruct.MaxGenFlow = 162.5; %flow in 1000 cfs
Dam(5).VariableStruct.MaxSpillFlow = 73.4; %flow in 1000 cfs
Dam(5).VariableStruct.MaxHead = 48.39; 
Dam(5).VariableStruct.MinHead = 30; %guess minimum height allowable

%% Double check max flow and head height. Max gen capacity is 121% of this flow & head
Dam(6).Type = 'Hydro Storage';
Dam(6).Name = 'Wanapum'; 
Dam(6).Source = 'Water';
Dam(6).Output = [];
Dam(6).Size = 163; %This would be storage capacity in kilo-acre-feet
Dam(6).Enabled = 1;
Dam(6).VariableStruct.MaxGenCapacity =  1517700; %power in kW
Dam(6).VariableStruct.RampUp =  924; %power change in kW/hr
Dam(6).VariableStruct.RampDown = 804; 
Dam(6).VariableStruct.MaxGenFlow = 173.7; %flow in 1000 cfs
Dam(6).VariableStruct.MaxSpillFlow = 97.6; %flow in 1000 cfs
Dam(6).VariableStruct.MaxHead = 85.03; 
Dam(6).VariableStruct.MinHead = 50; %guess minimum height allowable

%% Double check max flow and head height. Max gen capacity is 130% of this flow & head
Dam(7).Type = 'Hydro Storage';
Dam(7).Name = 'Priest Rapids';
Dam(7).Source = 'Water';
Dam(7).Output = [];
Dam(7).Size = 44; %This would be storage capacity in kilo-acre-feet
Dam(7).Enabled = 1;
Dam(7).VariableStruct.MaxGenCapacity = 1553600; %power in kW
Dam(7).VariableStruct.RampUp = 913; %power change in kW/hr
Dam(7).VariableStruct.RampDown = 809; 
Dam(7).VariableStruct.MaxGenFlow = 162.8; %flow in 1000 cfs
Dam(7).VariableStruct.MaxSpillFlow = 134.1; %flow in 1000 cfs
Dam(7).VariableStruct.MaxHead = 86.14; 
Dam(7).VariableStruct.MinHead = 40; %guess minimum height allowable

Dam(8).Type = 'Hydro Storage';
Dam(8).Name =  'McNary'; 
Dam(8).Source = 'Water';
Dam(8).Output = [];
Dam(8).Size = 1345; %This would be storage capacity in kilo-acre-feet
Dam(8).Enabled = 1;
Dam(8).VariableStruct.MaxGenCapacity = 1106540; %power in kW
Dam(8).VariableStruct.RampUp = 238; %power change in kW/hr
Dam(8).VariableStruct.RampDown = 378; 
Dam(8).VariableStruct.MaxGenFlow = 229.9; %flow in 1000 cfs
Dam(8).VariableStruct.MaxSpillFlow = 176.8; %flow in 1000 cfs
Dam(8).VariableStruct.MaxHead = 76.39; 
Dam(8).VariableStruct.MinHead = 40; %guess minimum height allowable

Dam(9).Type = 'Hydro Storage';
Dam(9).Name = 'John Day'; 
Dam(9).Source = 'Water';
Dam(9).Output = [];
Dam(9).Size = 2530; %This would be storage capacity in kilo-acre-feet
Dam(9).Enabled = 1;
Dam(9).VariableStruct.MaxGenCapacity = 1932770; %power in kW
Dam(9).VariableStruct.RampUp = 221; %power change in kW/hr
Dam(9).VariableStruct.RampDown = 223; 
Dam(9).VariableStruct.MaxGenFlow = 260.5; %flow in 1000 cfs
Dam(9).VariableStruct.MaxSpillFlow = 119; %flow in 1000 cfs
Dam(9).VariableStruct.MaxHead = 106.42; 
Dam(9).VariableStruct.MinHead = 60; %guess minimum height allowable

Dam(10).Type = 'Hydro Storage';
Dam(10).Name = 'The Dalles'; 
Dam(10).Source = 'Water';
Dam(10).Output = [];
Dam(10).Size = 330; %This would be storage capacity in kilo-acre-feet
Dam(10).Enabled = 1;
Dam(10).VariableStruct.MaxGenCapacity = 2160000; %power in kW
Dam(10).VariableStruct.RampUp = 227; %power change in kW/hr
Dam(10).VariableStruct.RampDown = 228; 
Dam(10).VariableStruct.MaxGenFlow = 263.8; %flow in 1000 cfs
Dam(10).VariableStruct.MaxSpillFlow = 128.1; %flow in 1000 cfs
Dam(10).VariableStruct.MaxHead = 85.79; 
Dam(10).VariableStruct.MinHead = 50; %guess minimum height allowable

Dam(11).Type = 'Hydro Storage';
Dam(11).Name = 'Bonneville'; 
Dam(11).Source = 'Water';
Dam(11).Output = [];
Dam(11).Size = 537; %This would be storage capacity in kilo-acre-feet
Dam(11).Enabled = 1;
Dam(11).VariableStruct.MaxGenCapacity = 1227000; %power in kW
Dam(11).VariableStruct.RampUp = 319; %power change in kW/hr
Dam(11).VariableStruct.RampDown = 370;  
Dam(11).VariableStruct.MaxGenFlow = 268.2; %flow in 1000 cfs
Dam(11).VariableStruct.MaxSpillFlow = 164.5; %flow in 1000 cfs
Dam(11).VariableStruct.MaxHead = 69.1; %ft
Dam(11).VariableStruct.MinHead = 40; %guess minimum height allowable

%% currently assuming 75% efficient to calculate max flow
Dam(12).Type = 'Hydro Storage';
Dam(12).Name = 'Brownlee'; 
Dam(12).Source = 'Water';
Dam(12).Output = [];
Dam(12).Size = 1426.7; %This would be storage capacity in kilo-acre-feet
Dam(12).Enabled = 1;
Dam(12).VariableStruct.MaxGenCapacity = 585400; %power in kW
Dam(12).VariableStruct.RampUp = 200; %power change in kW/hr
Dam(12).VariableStruct.RampDown = 200;
Dam(12).VariableStruct.MaxGenFlow = 22;
Dam(12).VariableStruct.MaxSpillFlow = 15;
Dam(12).VariableStruct.MaxHead = 420;
Dam(12).VariableStruct.MinHead = 325; %guess minimum height allowable

%% currently assuming 75% efficient to calculate max flow
Dam(13).Type = 'Hydro Storage';
Dam(13).Name = 'Oxbow'; 
Dam(13).Source = 'Water';
Dam(13).Output = [];
Dam(13).Size = 58.2; %This would be storage capacity in kilo-acre-feet
Dam(13).Enabled = 1;
Dam(13).VariableStruct.MaxGenCapacity = 190000; %power in kW
Dam(13).VariableStruct.RampUp = 200;%power change in kW/hr
Dam(13).VariableStruct.RampDown = 200;
Dam(13).VariableStruct.MaxGenFlow = 17; %flow in 1000 cfs
Dam(13).VariableStruct.MaxSpillFlow =15; %flow in 1000 cfs
Dam(13).VariableStruct.MaxHead = 175;
Dam(13).VariableStruct.MinHead = 100; %guess minimum height allowable

%% currently assuming 75% efficient to calculate max flow
Dam(14).Type = 'Hydro Storage';
Dam(14).Name = 'Hells Canyon';
Dam(14).Source = 'Water';
Dam(14).Output = [];
Dam(14).Size = 188; %This would be storage capacity in kilo-acre-feet
Dam(14).Enabled = 1;
Dam(14).VariableStruct.MaxGenCapacity = 391000; %power in kW
Dam(14).VariableStruct.RampUp = 200; %power change in kW/hr
Dam(14).VariableStruct.RampDown = 200;
Dam(14).VariableStruct.MaxGenFlow = 19; %flow in 1000 cfs
Dam(14).VariableStruct.MaxSpillFlow =15; %flow in 1000 cfs
Dam(14).VariableStruct.MaxHead = 330;
Dam(14).VariableStruct.MinHead = 180; %guess minimum height allowable
 
Dam(15).Type = 'Hydro Storage';
Dam(15).Name =  'Lower Granite'; 
Dam(15).Source = 'Water';
Dam(15).Output = [];
Dam(15).Size = 440; %This would be storage capacity in kilo-acre-feet
Dam(15).Enabled = 1;
Dam(15).VariableStruct.MaxGenCapacity = 627260; %power in kW
Dam(15).VariableStruct.RampUp = 215; %power change in kW/hr
Dam(15).VariableStruct.RampDown = 251; 
Dam(15).VariableStruct.MaxGenFlow = 89.6; %flow in 1000 cfs
Dam(15).VariableStruct.MaxSpillFlow = 46.5; %flow in 1000 cfs
Dam(15).VariableStruct.MaxHead = 102.5; 
Dam(15).VariableStruct.MinHead = 50; %guess minimum height allowable

Dam(16).Type = 'Hydro Storage';
Dam(16).Name = 'Little Goose'; 
Dam(16).Source = 'Water';
Dam(16).Output = [];
Dam(16).Size = 516; %This would be storage capacity in kilo-acre-feet
Dam(16).Enabled = 1;
Dam(16).VariableStruct.MaxGenCapacity = 642820; %power in kW
Dam(16).VariableStruct.RampUp = 373; %power change in kW/hr
Dam(16).VariableStruct.RampDown = 427; 
Dam(16).VariableStruct.MaxGenFlow = 91; %flow in 1000 cfs
Dam(16).VariableStruct.MaxSpillFlow = 38.1; %flow in 1000 cfs
Dam(16).VariableStruct.MaxHead = 98.96; 
Dam(16).VariableStruct.MinHead = 50; %guess minimum height allowable

Dam(17).Type = 'Hydro Storage';
Dam(17).Name = 'Lower Monumental';
Dam(17).Source = 'Water';
Dam(17).Output = [];
Dam(17).Size = 432; %This would be storage capacity in kilo-acre-feet
Dam(17).Enabled = 1;
Dam(17).VariableStruct.MaxGenCapacity = 699070; %power in kW
Dam(17).VariableStruct.RampUp = 274; %power change in kW/hr
Dam(17).VariableStruct.RampDown = 273; 
Dam(17).VariableStruct.MaxGenFlow = 95.3; %flow in 1000 cfs
Dam(17).VariableStruct.MaxSpillFlow = 47.1; %flow in 1000 cfs
Dam(17).VariableStruct.MaxHead = 441.41; 
Dam(17).VariableStruct.MinHead = 300; %guess minimum height allowable

Dam(18).Type = 'Hydro Storage';
Dam(18).Name =  'Ice Harbor';
Dam(18).Source = 'Water';
Dam(18).Output = [];
Dam(18).Size = 249; %This would be storage capacity in kilo-acre-feet
Dam(18).Enabled = 1;
Dam(18).VariableStruct.MaxGenCapacity = 581060; %power in kW
Dam(18).VariableStruct.RampUp = 609; %power change in kW/hr
Dam(18).VariableStruct.RampDown = 510;
Dam(18).VariableStruct.MaxGenFlow = 82.9; %flow in 1000 cfs
Dam(18).VariableStruct.MaxSpillFlow = 86.3; %flow in 1000 cfs
Dam(18).VariableStruct.MaxHead = 102.32;
Dam(18).VariableStruct.MinHead = 50; %guess minimum height allowable

DamNames = cell(length(Dam),1);
for i = 1:1:length(Dam)
    DamNames(i) = {Dam(i).Name};
end

Plant.Generator = [];
Plant.Network = [];
for i = 1:1:length(NodeNames);%number of dams in network
    Plant.Network(i).name = NodeNames{i}; %names of each node are name of closest city
    Plant.Network(i).Equipment = strcat('Hydro Storage.',Equipment(i)); %equiptment Hydro.Name_of_Dam
    
    Plant.Network(i).Electrical.connections = {};
    Plant.Network(i).Electrical.Trans_Eff = [];
    Plant.Network(i).Electrical.Trans_Limit = [];
    
    Plant.Network(i).Hydro.connections = {};
    Plant.Network(i).Hydro.Time2Sea = Time2Sea(i);
    Plant.Network(i).Hydro.InstreamFlow = InstreamFlow(i); %Specify instream flow requirements at the node point in the river
    up = UpRiver{i};
    for j = 1:1:length(up)
        if ~isempty(up{j})
            Plant.Network(i).Electrical.connections(end+1) = up(j);
            Plant.Network(i).Electrical.Trans_Eff(end+1) = 1;
            Plant.Network(i).Electrical.Trans_Limit(end+1) = inf;
        end
    end
    if ~isempty(DownRiver{i})%at most 1 downriver connection (zero for last dam before sea)  
        Plant.Network(i).Electrical.connections(end+1) = DownRiver(i);
        Plant.Network(i).Electrical.Trans_Eff(end+1) = 1;
        Plant.Network(i).Electrical.Trans_Limit(end+1) = inf;
        
        Plant.Network(i).Hydro.connections(end+1) = DownRiver(i);
    end
    I = find(strcmp(Equipment{i},DamNames));
    if ~isempty(I)
        Hydro_Eff(i) = Dam(I).VariableStruct.MaxGenCapacity/(Dam(I).VariableStruct.MaxGenFlow*Dam(I).VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW; 
        Dam(I).VariableStruct.RampUp = 1/5*Dam(I).VariableStruct.MaxGenCapacity;
        Dam(I).VariableStruct.RampDown = 1/5*Dam(I).VariableStruct.MaxGenCapacity;
        if i ==1
            Plant.Generator = Dam(I);
        else Plant.Generator(i) = Dam(I);
        end
    else
        %%put in irrigation district here
        I = find(strcmp(Equipment{i},IrrigationNames));
    end
end
Plant.Network(1).Electrical.Load = 1; %put all the load at the first node
% calculateHistoricalFit %might need to update this to do extra hydro specific fitting of data

%% Temporary filling in missing data
Plant.Data.Hydro.Timestamp = Plant.Data.Timestamp;
% Plant.Data.Hydro.Equipment = Equipment;
Plant.Data.Hydro.Nodes = NodeNames;
Plant.Data.Hydro.SourceSink(:,15:18) = Plant.Data.Hydro.SourceSink(:,13:16);
Plant.Data.Hydro.SourceSink(:,11:14) = 0;
for i = 1:1:11
    for k = 1:1:length(Plant.Data.Hydro.PowerGen(:,i))
        if isnan(Plant.Data.Hydro.PowerGen(k,i))
            Plant.Data.Hydro.PowerGen(k,i) = Plant.Data.Hydro.PowerGen(k-1,i);
        end
        if isnan(Plant.Data.Hydro.SourceSink(k,i))
            Plant.Data.Hydro.SourceSink(k,i) = Plant.Data.Hydro.SourceSink(k-1,i);
        end
        if isnan(Plant.Data.Hydro.InFlow(k,i))
            Plant.Data.Hydro.InFlow(k,i) = Plant.Data.Hydro.InFlow(k-1,i);
        end
    end
    HydroPerc(:,i) = Plant.Data.Hydro.PowerGen(:,i)/Plant.Generator(i).VariableStruct.MaxGenCapacity*100;
    HydroPercMean(:,i) = HydroPerc(:,i)/mean(HydroPerc(:,i));
end
for i = 1:1:4
    for k = 1:1:length(Plant.Data.Hydro.PowerGen(:,i+14))
        if isnan(Plant.Data.Hydro.PowerGen(k,i+14))
            Plant.Data.Hydro.PowerGen(k,i+14) = Plant.Data.Hydro.PowerGen(k-1,i+14);
        end
        if isnan(Plant.Data.Hydro.SourceSink(k,i+14))
            Plant.Data.Hydro.SourceSink(k,i+14) = Plant.Data.Hydro.SourceSink(k-1,i+14);
        end
        if isnan(Plant.Data.Hydro.InFlow(k,i+14))
            Plant.Data.Hydro.InFlow(k,i+14) = Plant.Data.Hydro.InFlow(k-1,i+14);
        end
    end
    HydroPerc2(:,i) = Plant.Data.Hydro.PowerGen(:,i+14)/Plant.Generator(i+14).VariableStruct.MaxGenCapacity*100;
    HydroPerc2Mean(:,i) = HydroPerc2(:,i)/mean(HydroPerc2(:,i));
end

for i = 1:1:3
    Plant.Data.Hydro.PowerGen(:,i+11) = round(mean(HydroPerc2,2)/100*Plant.Generator(i+11).VariableStruct.MaxGenCapacity/10)*10;
    Eff = Plant.Generator(i+11).VariableStruct.MaxGenCapacity/(Plant.Generator(i+11).VariableStruct.MaxGenFlow*Plant.Generator(i+11).VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW
    Plant.Data.Hydro.PowerFlow(:,i+11) = Plant.Data.Hydro.PowerGen(:,i+11)/(Eff*Plant.Generator(i+11).VariableStruct.MaxHead*84.674);%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 87.674 kJ/ (1000ft^3*ft)
    Plant.Data.Hydro.OutFlow(:,i+11) = Plant.Data.Hydro.PowerFlow(:,i+11);
    if i ==1
        Plant.Data.Hydro.SourceSink(:,12) = Plant.Data.Hydro.SourceSink(:,15) - 3*mean(Plant.Data.Hydro.SourceSink(:,16:18),2);
        Plant.Data.Hydro.InFlow(:,12) = Plant.Data.Hydro.SourceSink(:,12);
    end
    Plant.Data.Hydro.SourceSink(:,i+12) = mean(Plant.Data.Hydro.SourceSink(:,16:18),2);
    if i<3
        Plant.Data.Hydro.InFlow(:,i+12) = Plant.Data.Hydro.OutFlow(:,i+11);
    end
end
outputE = zeros(18,1);
for i = 1:1:18
    Eff = Plant.Generator(i).VariableStruct.MaxGenCapacity/(Plant.Generator(i).VariableStruct.MaxGenFlow*Plant.Generator(i).VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW
    outputE(i) = Eff*Plant.Generator(i).VariableStruct.MaxHead*84.674;%Power (kW) = efficiency(%) * Flow (1000 ft^3/s) * Head (ft) * 87.674 kJ/ (1000ft^3*ft)
end
Plant.Data.Demand.E = Plant.Data.Hydro.PowerFlow*outputE;