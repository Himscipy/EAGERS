global Plant
Plant = [];
Plant.Name = 'ColumbiaBasin';
%%Columbia river: 'GC';'CJ';'W';'RR';'RI';'WNPM';'PR';'MN';'JD';'TD';'BNVL';
%%Snake River: 'BRWNL';'OX';'HC';'LWRGRNT';'LTLGSE';'LWRMONUM';'ICEHRBR';
%%Excluding: BRWNL;OX;HC; We do not currently have enough Data

%%rows 1->366 (hourly): October 2nd, 2015 -> October 1st, 2016; 
Plant.Data.Timestamp(:,1)= linspace(datenum([2015 10 2 1 0 0]),datenum([2016 10 2 0 0 0]),8784).';
Plant.Data.Temperature = 20*ones(length(Plant.Data.Timestamp),1);
Plant.Data.Holidays = [];
Plant.Data.HistProf = [];

% %% Columns of SourceSink, SpillFlow, InFlow and Outflow must correspond to the node in the order of NodeNames
% %%columns 1->18 (kcfs or KW): Grand Coulee, Chief Joseph, Wells, Rocky Ridge, Rock Island, Wanapum, PriestRiver, McNary, John Day, The Dalles, Bonneville, Brownlee, Oxbow, Hells Canyon, Lower Granite, Little Goose, Lower Monumental, Ice Harbor;
% %SourceSink = Mass balance; Sinks and Sources in river segments; Negative values = sinks; Positive values = sources;
% load('SourcesandSinks'); load('Powerflow');load('Spillflow');load('Outflow');load('Inflow');load('PowerGen');
% Plant.Data.Hydro.SourceSink = SourcesandSinks; 
% %Powerflow = Flow used to produce power
% Plant.Data.Hydro.PowerFlow = Powerflow;
% %Spillflow = Flow that is spilled not used for power
% Plant.Data.Hydro.SpillFlow = Spillflow;
% %Outflow = Full amount of discharge flow
% Plant.Data.Hydro.OutFlow = Outflow;
% %columns (kilo-cfs): inflow(GC), infl(CJ)-Disch(GC), infl(W)-Disch(CJ), infl(RR)-Disch(W), infl(RI)-Disch(RR), infl(WNP)-Disch(RI), infl(PR)-Disch(WNP), infl(MN)-Disch(IH)-Disch(PR), infl(JD)-Disch(MN), infl(TD)-Disch(JD), infl(BNVL)-Disch(TD), inflow(BRW), infl(OX)-Disch(BRW), infl(HC)-Disch(OX), infl(LGRT)-Disch(HC), infl(LGS)-Disch(LGRT), infl(LMNT)-Disch(LG), infl(IH)-Disch(LMNT);
% Plant.Data.Hydro.InFlow = Inflow;
% %PowerGen = Generation (kW) every hour
% Plant.Data.Hydro.PowerGen = PowerGen;
% PowerGen(isnan(PowerGen)) = 0;
% Plant.Data.Demand = [];
% %Still need 
% %            storage data

%% Network Setup
%Brownlee, Oxbow, and Hells Canyon Dams ommitted

NodeNames = {'Plant1';'Plant2';'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa'};%'Plant3';'Plant4';'Plant5';'Plant6';'Plant7';'Plant8';...
%              'Plant9';'Plant10';'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa';'Wenatchee, Wa';...
%              'South Wenatchee, Wa';'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';...
%              'Rufus, Or';'The Dalles, Or';'Bonneville, Or';'Almota, Wa'; 'Starbuck, Wa';...
%              'Kahlotus, Wa';'Pasco, Wa';}; %15 nodes in the network + 11 Electric nodes 
         
DownRiver = {{'Plant2'};{''};'Bridgeport, Wa';'Azwell, Wa';''};
%                 'Plant2';'Plant3';'Plant4'};{'Plant3'};{'Plant4'};...
%              {'Plant5';'Plant7';'Plant8';'Plant9'};{'Plant6';'Plant9'};...
%              {'Plant7';'Plant9'};{'Plant8'};{''};{'Plant10'};{''};'Bridgeport, Wa';...
%              'Azwell, Wa';'Wenatchee, Wa';'South Wenatchee, Wa';...
%              'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';'Rufus, Or';...
%              'The Dalles, Or';'Bonneville, Or'; '';'Starbuck, Wa';...
%              'Kahlotus, Wa';'Pasco, Wa';''}; %Water connections; Umatilla, Or
           
UpRiver =  {{''};{'Plant1'}};%{'Plant1';'Plant2'};{'Plant1';'Plant3'};...
%             {'Plant4'};{'Plant5'};{'Plant4';'Plant6'};{'Plant4';'Plant7'};...
%             {'Plant4';'Plant5';'Plant6'};{'Plant9'}}; %Electrical network connections

Equipment = {'Plant1';'Plant2';'Grand Coulee';'Chief Joseph';'Wells'};%'Plant3';'Plant4';'Plant5';...
%              'Plant6';'Plant7';'Plant8';'Plant9';'Plant10';...
%              'Grand Coulee';'Chief Joseph';'Wells';'Rocky Reach';'Rock Island';...
%              'Wanapum';'Priest Rapids';'McNary';'John Day';'The Dalles';...
%              'Bonneville';'Lower Granite';'Little Goose';'Lower Monumental';...
%              'Ice Harbor'}; %Equipment at each node

HydroConnection = {'Plant1';'Plant2';'Plant2'};%'Plant3';'Plant3';'Plant4';'Plant4';'Plant6';...
%                 'Plant7';'Plant7';'Plant8';'Plant10';'Plant9';'Plant9';'Plant5'};             %BRWN/OX/HC: {'Cambridge, Id';'Oxbow, Or';'Hells Canyon, Or'}
         
InstreamFlow = [ 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;];
Time2Sea = [80.6;73.7;69.7;64.0;61.3;56.2;53.7;39.5;29.1;25.9;19.7;58.5;53.4;49.5;45.2];
RiverMile = [596.6;545.1;515.8;473.7;453.4;415.8;397.1;292;215.6;191.5;146.1;107.5;70.3;41.6;9.7;]; %was used to calculate the Time2Sea

%optimization options
Plant.optimoptions.Interval = 31;
Plant.optimoptions.Horizon = 24*7;
Plant.optimoptions.Resolution = 4;
Plant.optimoptions.Topt = 600;
Plant.optimoptions.Tmpc = 60;
Plant.optimoptions.nsSmooth = 0;
Plant.optimoptions.scaletime = 1;
Plant.optimoptions.fastsimulation = 1;
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
% Dam(1).VariableStruct.MaxGenCapacity = 4414000; %power in kW
% % Dam(1).VariableStruct.Capacityfactor = 0.36; 
% Dam(1).VariableStruct.RampUp = 189; %power change in kW/hr
% Dam(1).VariableStruct.RampDown = 233;
% Dam(1).VariableStruct.MaxGenFlow = 202.5; %flow in 1000 cfs
% Dam(1).VariableStruct.MaxSpillFlow = 0.2; %flow in 1000 cfs
% Dam(1).VariableStruct.MaxHead = 334.2;
% Dam(1).VariableStruct.MinHead = 180; %guess minimum height allowable
Dam(1).VariableStruct.MaxGenCapacity = 6079000; %power in kW
Dam(1).VariableStruct.RampUp = 3775000; %power change in kW/hr
Dam(1).VariableStruct.RampDown = 3938000;
Dam(1).VariableStruct.MaxGenFlow = 242.7; %flow in 1000 cfs
Dam(1).VariableStruct.MaxSpillFlow = 105.1; %flow in 1000 cfs
Dam(1).VariableStruct.MaxHead = 337;
Dam(1).VariableStruct.MinHead = 170; %guess minimum height allowable

Dam(2).Type = 'Hydro Storage';
Dam(2).Name = 'Chief Joseph'; 
Dam(2).Source = 'Water';
Dam(2).Output = [];
Dam(2).Size = 516; %This would be storage capacity in kilo-acre-feet
Dam(2).Enabled = 1;
% Dam(2).VariableStruct.MaxGenCapacity = 2180000; %power in kW
% Dam(2).VariableStruct.RampUp =532; %power change in kW/hr
% Dam(2).VariableStruct.RampDown =586; 
% Dam(2).VariableStruct.MaxGenFlow = 179.3; %flow in 1000 cfs
% Dam(2).VariableStruct.MaxSpillFlow = 31.3; %flow in 1000 cfs
% Dam(2).VariableStruct.MaxHead = 180.8; 
% Dam(2).VariableStruct.MinHead = 100; %guess minimum height allowable
Dam(2).VariableStruct.MaxGenCapacity = 2620000; %power in kW
Dam(2).VariableStruct.RampUp =1631000; %power change in kW/hr
Dam(2).VariableStruct.RampDown =1683000; 
Dam(2).VariableStruct.MaxGenFlow = 212.8; %flow in 1000 cfs
Dam(2).VariableStruct.MaxSpillFlow = 196.5; %flow in 1000 cfs
Dam(2).VariableStruct.MaxHead = 199.8; 
Dam(2).VariableStruct.MinHead = 159; %guess minimum height allowable

Dam(3).Type = 'Hydro Storage';
Dam(3).Name = 'Wells';
Dam(3).Source = 'Water';
Dam(3).Output = [];
Dam(3).Size = 331.2; %This would be storage capacity in kilo-acre-feet
% Dam(3).Size = 98; %This would be storage capacity in kilo-acre-feet
Dam(3).Enabled = 1;
% Dam(3).VariableStruct.MaxGenCapacity = 851400; %power in kW
% Dam(3).VariableStruct.RampUp =533; %power change in kW/hr
% Dam(3).VariableStruct.RampDown =608; 
% Dam(3).VariableStruct.MaxGenFlow = 175.3; %flow in 1000 cfs
% Dam(3).VariableStruct.MaxSpillFlow = 57.96; %flow in 1000 cfs
% Dam(3).VariableStruct.MaxHead = 75.19; 
% Dam(3).VariableStruct.MinHead = 40; %guess minimum height allowable
Dam(3).VariableStruct.MaxGenCapacity = 851400; %power in kW
Dam(3).VariableStruct.RampUp =599318; %power change in kW/hr
Dam(3).VariableStruct.RampDown =613318; 
Dam(3).VariableStruct.MaxGenFlow = 219.2; %flow in 1000 cfs
Dam(3).VariableStruct.MaxSpillFlow = 185.5; %flow in 1000 cfs
Dam(3).VariableStruct.MaxHead = 76.65; 
Dam(3).VariableStruct.MinHead = 46.91; %guess minimum height allowable
% 
% Dam(4).Type = 'Hydro Storage';
% Dam(4).Name =  'Rocky Reach'; 
% Dam(4).Source = 'Water';
% Dam(4).Output = [];
% Dam(4).Size = 382; %This would be storage capacity in kilo-acre-feet
% % Dam(4).Size = 37; %This would be storage capacity in kilo-acre-feet
% Dam(4).Enabled = 1;
% % Dam(4).VariableStruct.MaxGenCapacity = 1012000; %power in kW
% % Dam(4).VariableStruct.RampUp =396; %power change in kW/hr
% % Dam(4).VariableStruct.RampDown =476; 
% % Dam(4).VariableStruct.MaxGenFlow = 159.8; %flow in 1000 cfs
% % Dam(4).VariableStruct.MaxSpillFlow = 52.99; %flow in 1000 cfs
% % Dam(4).VariableStruct.MaxHead =  95.19; 
% % Dam(4).VariableStruct.MinHead = 60; %guess minimum height allowable
% Dam(4).VariableStruct.MaxGenCapacity = 1299600; %power in kW
% Dam(4).VariableStruct.RampUp =705000; %power change in kW/hr
% Dam(4).VariableStruct.RampDown =632000; 
% Dam(4).VariableStruct.MaxGenFlow = 365.4; %flow in 1000 cfs
% Dam(4).VariableStruct.MaxSpillFlow = 200.1; %flow in 1000 cfs
% Dam(4).VariableStruct.MaxHead =  97.81; 
% Dam(4).VariableStruct.MinHead = 73.56; %guess minimum height allowable
%  
% Dam(5).Type = 'Hydro Storage';
% Dam(5).Name = 'Rock Island'; 
% Dam(5).Source = 'Water';
% Dam(5).Output = [];
% Dam(5).Size = 131; %This would be storage capacity in kilo-acre-feet
% Dam(5).Enabled = 1;
% Dam(5).VariableStruct.MaxGenCapacity = 623700; %power in kW
% % Dam(5).VariableStruct.RampUp = 391; %power change in kW/hr
% % Dam(5).VariableStruct.RampDown = 478; 
% % Dam(5).VariableStruct.MaxGenFlow = 162.5; %flow in 1000 cfs
% % Dam(5).VariableStruct.MaxSpillFlow = 73.4; %flow in 1000 cfs
% % Dam(5).VariableStruct.MaxHead = 48.39; 
% % Dam(5).VariableStruct.MinHead = 30; %guess minimum height allowable
% Dam(5).VariableStruct.RampUp = 342000; %power change in kW/hr
% Dam(5).VariableStruct.RampDown = 307000; 
% Dam(5).VariableStruct.MaxGenFlow = 235.6; %flow in 1000 cfs
% Dam(5).VariableStruct.MaxSpillFlow = 153.7; %flow in 1000 cfs
% Dam(5).VariableStruct.MaxHead = 50.76; 
% Dam(5).VariableStruct.MinHead = 32.21; %guess minimum height allowable
% 
% Dam(6).Type = 'Hydro Storage';
% Dam(6).Name = 'Wanapum'; 
% Dam(6).Source = 'Water';
% Dam(6).Output = [];
% Dam(6).Size = 796; %This would be storage capacity in kilo-acre-feet
% Dam(6).Enabled = 1;
% Dam(6).VariableStruct.MaxGenCapacity = 1040000; %power in kW
% % Dam(6).VariableStruct.RampUp =  924; %power change in kW/hr
% % Dam(6).VariableStruct.RampDown = 804; 
% % Dam(6).VariableStruct.MaxGenFlow = 173.7; %flow in 1000 cfs
% % Dam(6).VariableStruct.MaxSpillFlow = 97.6; %flow in 1000 cfs
% % Dam(6).VariableStruct.MaxHead = 85.03; 
% % Dam(6).VariableStruct.MinHead = 50; %guess minimum height allowable
% Dam(6).VariableStruct.RampUp =  444920; %power change in kW/hr
% Dam(6).VariableStruct.RampDown = 444910; 
% Dam(6).VariableStruct.MaxGenFlow = 192.7; %flow in 1000 cfs
% Dam(6).VariableStruct.MaxSpillFlow = 294.7; %flow in 1000 cfs
% Dam(6).VariableStruct.MaxHead = 91.29; 
% Dam(6).VariableStruct.MinHead = 44.17; %guess minimum height allowable
% 
% Dam(7).Type = 'Hydro Storage';
% Dam(7).Name = 'Priest Rapids';
% Dam(7).Source = 'Water';
% Dam(7).Output = [];
% Dam(7).Size = 237.1; %This would be storage capacity in kilo-acre-feet
% Dam(7).Enabled = 1;
% Dam(7).VariableStruct.MaxGenCapacity = 955600; %power in kW
% % Dam(7).VariableStruct.RampUp = 913; %power change in kW/hr
% % Dam(7).VariableStruct.RampDown = 809; 
% % Dam(7).VariableStruct.MaxGenFlow = 180.4; %flow in 1000 cfs
% % Dam(7).VariableStruct.MaxSpillFlow = 166.1; %flow in 1000 cfs
% % Dam(7).VariableStruct.MaxHead = 86.14; 
% % Dam(7).VariableStruct.MinHead = 40; %guess minimum height allowable
% Dam(7).VariableStruct.RampUp = 526350; %power change in kW/hr
% Dam(7).VariableStruct.RampDown = 513700; 
% Dam(7).VariableStruct.MaxGenFlow = 180.4; %flow in 1000 cfs
% Dam(7).VariableStruct.MaxSpillFlow = 306.1; %flow in 1000 cfs
% Dam(7).VariableStruct.MaxHead = 100.53; 
% Dam(7).VariableStruct.MinHead = 54.84;
% 
% Dam(8).Type = 'Hydro Storage';
% Dam(8).Name =  'McNary'; 
% Dam(8).Source = 'Water';
% Dam(8).Output = [];
% Dam(8).Size = 1350; %This would be storage capacity in kilo-acre-feet
% Dam(8).Enabled = 1;
% Dam(8).VariableStruct.MaxGenCapacity = 1127000; %power in kW
% % Dam(8).VariableStruct.RampUp = 238; %power change in kW/hr
% % Dam(8).VariableStruct.RampDown = 378; 
% % Dam(8).VariableStruct.MaxGenFlow = 229.9; %flow in 1000 cfs
% % Dam(8).VariableStruct.MaxSpillFlow = 176.8; %flow in 1000 cfs
% % Dam(8).VariableStruct.MaxHead = 76.39; 
% % Dam(8).VariableStruct.MinHead = 40; %guess minimum height allowable
% Dam(8).VariableStruct.RampUp = 355200; %power change in kW/hr
% Dam(8).VariableStruct.RampDown = 395330; 
% Dam(8).VariableStruct.MaxGenFlow = 231.1; %flow in 1000 cfs
% Dam(8).VariableStruct.MaxSpillFlow = 375.5; %flow in 1000 cfs
% Dam(8).VariableStruct.MaxHead = 78.14; 
% Dam(8).VariableStruct.MinHead = 63.9; %guess minimum height allowable
% 
% Dam(9).Type = 'Hydro Storage';
% Dam(9).Name = 'John Day'; 
% Dam(9).Source = 'Water';
% Dam(9).Output = [];
% Dam(9).Size = 2530; %This would be storage capacity in kilo-acre-feet
% Dam(9).Enabled = 1;
% % Dam(9).VariableStruct.MaxGenCapacity =1932770; %power in kW
% % Dam(9).VariableStruct.RampUp = 221; %power change in kW/hr
% % Dam(9).VariableStruct.RampDown = 223; 
% % Dam(9).VariableStruct.MaxGenFlow = 260.5; %flow in 1000 cfs
% % Dam(9).VariableStruct.MaxSpillFlow = 119; %flow in 1000 cfs
% % Dam(9).VariableStruct.MaxHead = 106.42; 
% % Dam(9).VariableStruct.MinHead = 60; %guess minimum height allowable
% Dam(9).VariableStruct.MaxGenCapacity = 2160000; %power in kW
% Dam(9).VariableStruct.RampUp = 1252000; %power change in kW/hr
% Dam(9).VariableStruct.RampDown = 1107000; 
% Dam(9).VariableStruct.MaxGenFlow = 319.9; %flow in 1000 cfs
% Dam(9).VariableStruct.MaxSpillFlow = 270.6; %flow in 1000 cfs
% Dam(9).VariableStruct.MaxHead = 112.9; 
% Dam(9).VariableStruct.MinHead = 92; %guess minimum height allowable
% 
% 
% Dam(10).Type = 'Hydro Storage';
% Dam(10).Name = 'The Dalles'; 
% Dam(10).Source = 'Water';
% Dam(10).Output = [];
% Dam(10).Size = 330; %This would be storage capacity in kilo-acre-feet
% Dam(10).Enabled = 1;
% % Dam(10).VariableStruct.MaxGenCapacity = 1878300; %power in kW
% % Dam(10).VariableStruct.RampUp = 227; %power change in kW/hr
% % Dam(10).VariableStruct.RampDown = 228; 
% % Dam(10).VariableStruct.MaxGenFlow = 334.8; %flow in 1000 cfs
% % Dam(10).VariableStruct.MaxSpillFlow = 310.2; %flow in 1000 cfs
% % Dam(10).VariableStruct.MaxHead = 85.79; 
% % Dam(10).VariableStruct.MinHead = 50; %guess minimum height allowable
% Dam(10).VariableStruct.MaxGenCapacity = 2160000; %power in kW
% Dam(10).VariableStruct.RampUp = 607000; %power change in kW/hr
% Dam(10).VariableStruct.RampDown = 657800; 
% Dam(10).VariableStruct.MaxGenFlow = 334.8; %flow in 1000 cfs
% Dam(10).VariableStruct.MaxSpillFlow = 336.2; %flow in 1000 cfs
% Dam(10).VariableStruct.MaxHead = 87.1; 
% Dam(10).VariableStruct.MinHead = 67.98; %guess minimum height allowable
% 
% Dam(11).Type = 'Hydro Storage';
% Dam(11).Name = 'Bonneville'; 
% Dam(11).Source = 'Water';
% Dam(11).Output = [];
% Dam(11).Size = 537; %This would be storage capacity in kilo-acre-feet
% Dam(11).Enabled = 1;
% Dam(11).VariableStruct.MaxGenCapacity = 1242000; %power in kW
% % Dam(11).VariableStruct.RampUp = 319; %power change in kW/hr
% % Dam(11).VariableStruct.RampDown = 370;  
% % Dam(11).VariableStruct.MaxGenFlow = 268.2; %flow in 1000 cfs
% % Dam(11).VariableStruct.MaxSpillFlow = 164.5; %flow in 1000 cfs
% % Dam(11).VariableStruct.MaxHead = 69.1; %ft
% % Dam(11).VariableStruct.MinHead = 40; %guess minimum height allowable
% Dam(11).VariableStruct.RampUp = 537000; %power change in kW/hr
% Dam(11).VariableStruct.RampDown = 491000;  
% Dam(11).VariableStruct.MaxGenFlow = 274.7; %flow in 1000 cfs
% Dam(11).VariableStruct.MaxSpillFlow = 317.8; %flow in 1000 cfs
% Dam(11).VariableStruct.MaxHead = 70.3; %ft
% Dam(11).VariableStruct.MinHead = 37; %guess minimum height allowable
% 
% %%Not Enough Data to use these Dams in simulation
% % %% currently assuming 75% efficient to calculate max flow
% % Dam(12).Type = 'Hydro Storage';
% % Dam(12).Name = 'Brownlee'; 
% % Dam(12).Source = 'Water';
% % Dam(12).Output = [];
% % Dam(12).Size = 1426.7; %This would be storage capacity in kilo-acre-feet
% % Dam(12).Enabled = 1;
% % Dam(12).VariableStruct.MaxGenCapacity = 585400; %power in kW
% % Dam(12).VariableStruct.RampUp = 200; %power change in kW/hr
% % Dam(12).VariableStruct.RampDown = 200;
% % Dam(12).VariableStruct.MaxGenFlow = 22;
% % Dam(12).VariableStruct.MaxSpillFlow = 15;
% % Dam(12).VariableStruct.MaxHead = 420;
% % Dam(12).VariableStruct.MinHead = 325; %guess minimum height allowable
% % 
% % %% currently assuming 75% efficient to calculate max flow
% % Dam(13).Type = 'Hydro Storage';
% % Dam(13).Name = 'Oxbow'; 
% % Dam(13).Source = 'Water';
% % Dam(13).Output = [];
% % Dam(13).Size = 58.2; %This would be storage capacity in kilo-acre-feet
% % Dam(13).Enabled = 1;
% % Dam(13).VariableStruct.MaxGenCapacity = 190000; %power in kW
% % Dam(13).VariableStruct.RampUp = 200;%power change in kW/hr
% % Dam(13).VariableStruct.RampDown = 200;
% % Dam(13).VariableStruct.MaxGenFlow = 17; %flow in 1000 cfs
% % Dam(13).VariableStruct.MaxSpillFlow =15; %flow in 1000 cfs
% % Dam(13).VariableStruct.MaxHead = 175;
% % Dam(13).VariableStruct.MinHead = 100; %guess minimum height allowable
% % 
% % %% currently assuming 75% efficient to calculate max flow
% % Dam(14).Type = 'Hydro Storage';
% % Dam(14).Name = 'Hells Canyon';
% % Dam(14).Source = 'Water';
% % Dam(14).Output = [];
% % Dam(14).Size = 188; %This would be storage capacity in kilo-acre-feet
% % Dam(14).Enabled = 1;
% % Dam(14).VariableStruct.MaxGenCapacity = 391000; %power in kW
% % Dam(14).VariableStruct.RampUp = 200; %power change in kW/hr
% % Dam(14).VariableStruct.RampDown = 200;
% % Dam(14).VariableStruct.MaxGenFlow = 19; %flow in 1000 cfs
% % Dam(14).VariableStruct.MaxSpillFlow =15; %flow in 1000 cfs
% % Dam(14).VariableStruct.MaxHead = 330;
% % Dam(14).VariableStruct.MinHead = 180; %guess minimum height allowable
% 
% Dam(12).Type = 'Hydro Storage';
% Dam(12).Name =  'Lower Granite'; 
% Dam(12).Source = 'Water';
% Dam(12).Output = [];
% Dam(12).Size = 440.2; %This would be storage capacity in kilo-acre-feet
% Dam(12).Enabled = 1;
% % Dam(12).VariableStruct.MaxGenCapacity = 627260; %power in kW
% % Dam(12).VariableStruct.RampUp = 215; %power change in kW/hr
% % Dam(12).VariableStruct.RampDown = 251; 
% % Dam(12).VariableStruct.MaxGenFlow = 89.6; %flow in 1000 cfs
% % Dam(12).VariableStruct.MaxSpillFlow = 46.5; %flow in 1000 cfs
% % Dam(12).VariableStruct.MaxHead = 102.5; 
% % Dam(12).VariableStruct.MinHead = 50; %guess minimum height allowable
% Dam(12).VariableStruct.MaxGenCapacity = 848060; %power in kW
% Dam(12).VariableStruct.RampUp = 336770; %power change in kW/hr
% Dam(12).VariableStruct.RampDown = 432000; 
% Dam(12).VariableStruct.MaxGenFlow = 122.9; %flow in 1000 cfs
% Dam(12).VariableStruct.MaxSpillFlow = 143.6; %flow in 1000 cfs
% Dam(12).VariableStruct.MaxHead = 104; 
% Dam(12).VariableStruct.MinHead = 93.15; %guess minimum height allowable
% 
% Dam(13).Type = 'Hydro Storage';
% Dam(13).Name = 'Little Goose'; 
% Dam(13).Source = 'Water';
% Dam(13).Output = [];
% Dam(13).Size = 516.3; %This would be storage capacity in kilo-acre-feet
% Dam(13).Enabled = 1;
% % Dam(13).VariableStruct.MaxGenCapacity = 642820; %power in kW
% % Dam(13).VariableStruct.RampUp = 373; %power change in kW/hr
% % Dam(13).VariableStruct.RampDown = 427; 
% % Dam(13).VariableStruct.MaxGenFlow = 91; %flow in 1000 cfs
% % Dam(13).VariableStruct.MaxSpillFlow = 38.1; %flow in 1000 cfs
% % Dam(13).VariableStruct.MaxHead = 98.96; 
% % Dam(13).VariableStruct.MinHead = 50; %guess minimum height allowable
% Dam(13).VariableStruct.MaxGenCapacity = 932000; %power in kW
% Dam(13).VariableStruct.RampUp = 517640; %power change in kW/hr
% Dam(13).VariableStruct.RampDown = 495170; 
% Dam(13).VariableStruct.MaxGenFlow = 122.4; %flow in 1000 cfs
% Dam(13).VariableStruct.MaxSpillFlow = 199.6; %flow in 1000 cfs
% Dam(13).VariableStruct.MaxHead = 100.3; 
% Dam(13).VariableStruct.MinHead = 90.32; %guess minimum height allowable
% 
% Dam(14).Type = 'Hydro Storage';
% Dam(14).Name = 'Lower Monumental';
% Dam(14).Source = 'Water';
% Dam(14).Output = [];
% Dam(14).Size = 432; %This would be storage capacity in kilo-acre-feet
% Dam(14).Enabled = 1;
% % Dam(14).VariableStruct.MaxGenCapacity = 699070; %power in kW
% % Dam(14).VariableStruct.RampUp = 274; %power change in kW/hr
% % Dam(14).VariableStruct.RampDown = 273; 
% % Dam(14).VariableStruct.MaxGenFlow = 95.3; %flow in 1000 cfs
% % Dam(14).VariableStruct.MaxSpillFlow = 47.1; %flow in 1000 cfs
% % Dam(14).VariableStruct.MaxHead = 103.6; 
% % Dam(14).VariableStruct.MinHead = 70.4; %guess minimum height allowable
% Dam(14).VariableStruct.MaxGenCapacity = 932000; %power in kW
% Dam(14).VariableStruct.RampUp = 536060; %power change in kW/hr
% Dam(14).VariableStruct.RampDown = 469440; 
% Dam(14).VariableStruct.MaxGenFlow = 117.3; %flow in 1000 cfs
% Dam(14).VariableStruct.MaxSpillFlow = 157.4; %flow in 1000 cfs
% Dam(14).VariableStruct.MaxHead = 103.54; 
% Dam(14).VariableStruct.MinHead = 89.72; %guess minimum height allowable
% 
% Dam(15).Type = 'Hydro Storage';
% Dam(15).Name =  'Ice Harbor';
% Dam(15).Source = 'Water';
% Dam(15).Output = [];
% Dam(15).Size = 249; %This would be storage capacity in kilo-acre-feet
% Dam(15).Enabled = 1;
% % Dam(15).VariableStruct.MaxGenCapacity = 581060; %power in kW
% % Dam(15).VariableStruct.RampUp = 609; %power change in kW/hr
% % Dam(15).VariableStruct.RampDown = 510;
% % Dam(15).VariableStruct.MaxGenFlow = 82.9; %flow in 1000 cfs
% % Dam(15).VariableStruct.MaxSpillFlow = 86.3; %flow in 1000 cfs
% % Dam(15).VariableStruct.MaxHead = 102.32;
% % Dam(15).VariableStruct.MinHead = 50; %guess minimum height allowable
% Dam(15).VariableStruct.MaxGenCapacity = 693000; %power in kW
% Dam(15).VariableStruct.RampUp = 401160; %power change in kW/hr
% Dam(15).VariableStruct.RampDown = 439460;
% Dam(15).VariableStruct.MaxGenFlow = 98.2; %flow in 1000 cfs
% Dam(15).VariableStruct.MaxSpillFlow = 175.1; %flow in 1000 cfs
% Dam(15).VariableStruct.MaxHead = 102.32;
% Dam(15).VariableStruct.MinHead = 39.3; %guess minimum height allowable

Electric(1).Type = 'Electric';
Electric(1).Name = 'Station1';
Electric(1).Source = 'Electric';
Electric(1).Output = [];
Electric(1).Size = [];
Electric(1).Enabled = 1;

Electric(2).Type = 'Electric';
Electric(2).Name = 'Station2';
Electric(2).Source = 'Electric';
Electric(2).Output = [];
Electric(2).Size = [];
Electric(2).Enabled = 1;

% Electric(3).Type = 'Electric';
% Electric(3).Name = 'Station3';
% Electric(3).Source = 'Electric';
% Electric(3).Output = [];
% Electric(3).Size = [];
% Electric(3).Enabled = 1;
% 
% Electric(4).Type = 'Electric';
% Electric(4).Name = 'Station4';
% Electric(4).Source = 'Electric';
% Electric(4).Output = [];
% Electric(4).Size = [];
% Electric(4).Enabled = 1;
% 
% Electric(5).Type = 'Electric';
% Electric(5).Name = 'Station5';
% Electric(5).Source = 'Electric';
% Electric(5).Output = [];
% Electric(5).Size = [];
% Electric(5).Enabled = 1;
% 
% Electric(6).Type = 'Electric';
% Electric(6).Name = 'Station6';
% Electric(6).Source = 'Electric';
% Electric(6).Output = [];
% Electric(6).Size = [];
% Electric(6).Enabled = 1;
% 
% Electric(7).Type = 'Electric';
% Electric(7).Name = 'Station7';
% Electric(7).Source = 'Electric';
% Electric(7).Output = [];
% Electric(7).Size = [];
% Electric(7).Enabled = 1;
% 
% Electric(8).Type = 'Electric';
% Electric(8).Name = 'Station8';
% Electric(8).Source = 'Electric';
% Electric(8).Output = [];
% Electric(8).Size = [];
% Electric(8).Enabled = 1;
% 
% Electric(9).Type = 'Electric';
% Electric(9).Name = 'Station9';
% Electric(9).Source = 'Electric';
% Electric(9).Output = [];
% Electric(9).Size = [];
% Electric(9).Enabled = 1;
% 
% Electric(10).Type = 'Electric';
% Electric(10).Name = 'Station10';
% Electric(10).Source = 'Electric';
% Electric(10).Output = [];
% Electric(10).Size = [];
% Electric(10).Enabled = 1;

%Node Contains: BRWN/OX/HC (NOT USED)
% Electric(11).Type = 'Electric';
% Electric(11).Name = 'Station11';
% Electric(11).Source = 'Electric';
% Electric(11).Output = [];
% Electric(11).Size = [];
% Electric(11).Enabled = 1;

DamNames = cell(length(Dam),1);
x = 0;
for i = 1:1:length(NodeNames)
    if isempty(strfind(NodeNames{i},'Plant'))
        DamNames(i-x) = {Dam(i-x).Name};
    else 
        x = x + 1;
    end
end


%% Loading Data
load('allHistHydroGenData_2007_2016'); %Historical Generation Data(hourly)(kW)
% load('SourcesandSinks_2007_2016'); %Source and Sink Data; UpstreamOutFlow(t-T)-DownstreamInFlow(t)
%%columns 1->18 (kcfs or KW): Grand Coulee, Chief Joseph, Wells, Rocky Ridge, Rock Island, Wanapum, PriestRiver, McNary, John Day, The Dalles, Bonneville, Brownlee, Oxbow, Hells Canyon, Lower Granite, Little Goose, Lower Monumental, Ice Harbor;
%SourceSink = Mass balance; Sinks and Sources in river segments; Negative values = sinks; Positive values = sources;

%These are corrected data 
load('SourceSinkHrly')
load('InFlowCrctd')
load('OutFlowNMD')
load('powerFlowNMD')

for i = 1:1:length(Dam) 
    if exist(strcat(DamNames{i},'_','2007','_','2016','.mat'),'file')
        name = strcat(DamNames{i},'_','2007','_','2016');
        R = load(name);
        ix = fieldnames(R);
        Plant.Data.Hydro.Dams(i,1) = { R.(ix{1}) };
    else Plant.Data.Hydro.Dams{i,1} = {};
    end
end 

%10/2/15 1:00am thorugh 10/2/16 0:00am
firstday = datenum([2015 10 2 1 0 0]); %42279 in matlab
lastday = datenum([2016 10 2 1 0 0]);  %42645 in matlab
adjDate = datenum([1899 12 30 1 0 0]);
StartDate = firstday-adjDate;
EndDate = lastday-adjDate;
xi = [];
xf = [];
%Matlab treats a new day as midnight, while my data treats a new day as 1am
%Bringing data into alignment with Matlab
for i = 1:1:length(Dam)
    for j = 1:1:length(Plant.Data.Hydro.Dams{i})
        if Plant.Data.Hydro.Dams{i}(j,2) == 2400
            Plant.Data.Hydro.Dams{i}(j,2) =0;
            Plant.Data.Hydro.Dams{i}(j,1) = Plant.Data.Hydro.Dams{i}(j,1) + 1;
        end
        if Plant.Data.Hydro.Dams{i}(j,1) == StartDate && Plant.Data.Hydro.Dams{i}(j,2) == 100 && isempty(xi)
            xi = j;
        elseif Plant.Data.Hydro.Dams{i}(j,1) == EndDate  && isempty(xf)
            xf = j;
        end
    end
end 

for i = 1:1:length(Dam)
    Plant.Data.Hydro.SpillFlow(1:(xf-xi)+1,i) = Plant.Data.Hydro.Dams{i}(xi:xf,5);
    Plant.Data.Hydro.PowerGen(1:(xf-xi)+1,i) = Plant.Data.Hydro.Dams{i}(xi:xf,6);
    Plant.Data.Hydro.Head(1:(xf-xi)+1,i) = Plant.Data.Hydro.Dams{i}(xi:xf,7);
    if i == 1
        Plant.Data.Hydro.OutFlow(1:(xf-xi)+1,:) = OutFlow(xi:xf,1:length(Dam));
        Plant.Data.Hydro.InFlow(1:(xf-xi)+1,:) = crctdInFlow(xi:xf,1:length(Dam));
        Plant.Data.Hydro.PowerFlow(1:(xf-xi)+1,:) = PowerFlow(xi:xf,1:length(Dam));
        Plant.Data.Hydro.SourceSink(1:(xf-xi)+1,:) = SourceSinkHrly(xi:xf,1:length(Dam));
%         Plant.Data.Hydro.OutFlow(1:(xf-xi)+1,:) = OutFlow(xi:xf,:);
%         Plant.Data.Hydro.InFlow(1:(xf-xi)+1,:) = crctdInFlow(xi:xf,:);
%         Plant.Data.Hydro.PowerFlow(1:(xf-xi)+1,:) = PowerFlow(xi:xf,:);
%         Plant.Data.Hydro.SourceSink(1:(xf-xi)+1,:) = SourceSinkHrly(xi:xf,:);
    end 
    if i == 8
        Plant.Data.Hydro.SourceSink(1:(xf-xi)+1,i) = Plant.Data.Hydro.SourceSink(1:(xf-xi)+1,i)+Plant.Data.Hydro.OutFlow(1:(xf-xi)+1,15);
    end  
end 

%% Network Creation  

Plant.Generator = [];
Plant.Network = [];
Plant.Data.Hydro.Nodes = {};
for i = 1:1:length(NodeNames)%number of dams in network
    Plant.Network(i).name = NodeNames{i}; %names of each node are name of closest city
    Down = DownRiver{i};
    I = find(strcmp(Equipment{i},DamNames));
    Plant.Network(i).Equipment = {};
    if ~isempty(I)
        Plant.Data.Hydro.Nodes(end+1,1) = NodeNames(i);
        Plant.Network(i).Equipment = strcat('Hydro Storage.',Equipment(i)); %equiptment Hydro.Name_of_Dam   
        Plant.Network(i).Hydro.connections = {};
        Plant.Network(i).Hydro.Time2Sea = Time2Sea(I);
        Plant.Network(i).Hydro.InstreamFlow = InstreamFlow(I); %Specify instream flow requirements at the node point in the river
        if ~isempty(DownRiver{i})%at most 1 downriver connection (zero for last dam before sea)  
            Plant.Network(i).Hydro.connections(end+1) = DownRiver(i);
        end
        if ~isempty(I)
            Hydro_Eff(i) = Dam(I).VariableStruct.MaxGenCapacity/(Dam(I).VariableStruct.MaxGenFlow*Dam(I).VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW; 
            Dam(I).VariableStruct.RampUp = 1/5*Dam(I).VariableStruct.MaxGenCapacity;
            Dam(I).VariableStruct.RampDown = 1/5*Dam(I).VariableStruct.MaxGenCapacity;
            if I ==1
                Plant.Generator = Dam(I);
            else Plant.Generator(I) = Dam(I);
            end
        else
            %%put in irrigation district here
            I = find(strcmp(Equipment{i},IrrigationNames));
        end
        %Hydro => Plant power connections
        hydroConnect = HydroConnection{I};
        if ~isempty(hydroConnect)
            for k = 1:1:x
                if strcmp(hydroConnect,Plant.Network(k).name)
                    Plant.Network(i).Electrical.connections = {};
                    Plant.Network(i).Electrical.Trans_Eff = [];
                    Plant.Network(i).Electrical.Trans_Limit = [];
                    Plant.Network(i).Electrical.connections(end+1) = {Plant.Network(k).name};
                    Plant.Network(i).Electrical.Trans_Eff(end+1) = 1;
                    Plant.Network(i).Electrical.Trans_Limit(end+1) = inf;

                    %Adding other direction (dir = reverse) for electrical to dam
                    Plant.Network(k).Electrical.connections(end+1) = {Plant.Network(i).name};
                    Plant.Network(k).Electrical.Trans_Eff(end+1) = 1;
                    Plant.Network(k).Electrical.Trans_Limit(end+1) = inf;
                end  
            end 
        end 
    else
        up = UpRiver{i};
        Plant.Network(i).Electrical.connections = {};
        Plant.Network(i).Electrical.Trans_Eff = [];
        Plant.Network(i).Electrical.Trans_Limit = [];
%         Plant.Network(i).Equipment = strcat('Electric.',Equipment(i));
        for j = 1:1:length(Down)
            if ~isempty(Down{j})
                Plant.Network(i).Electrical.connections(end+1) = (Down(j));
                Plant.Network(i).Electrical.Trans_Eff(end+1) = 1-(i*.01);
                Plant.Network(i).Electrical.Trans_Limit(end+1) = inf;
            end  
        end
        for j = 1:1:length(up)
            if ~isempty(up{j})
                Plant.Network(i).Electrical.connections(end+1) = (up(j)); %upconnections electric
                Plant.Network(i).Electrical.Trans_Eff(end+1) = 1-(i*.01);
                Plant.Network(i).Electrical.Trans_Limit(end+1) = inf;
            end
        end 
    end 
end
Plant.Network(5).Electrical.Load = 1; %put all the load at the first node
% calculateHistoricalFit %might need to update this to do extra hydro specific fitting of data

%% Temporary filling in missing data
Plant.Data.Hydro.Timestamp = Plant.Data.Timestamp;
% Plant.Data.Hydro.Equipment = Equipment;
% Plant.Data.Hydro.Nodes = NodeNames{x+1:x+length(Dam)};
%%This was for: sourceBrwn,Ox-Brwn,HC-Ox,LG-HC
% Plant.Data.Hydro.SourceSink(:,15:18) = Plant.Data.Hydro.SourceSink(:,13:16);
% Plant.Data.Hydro.SourceSink(:,11:14) = 0;

%For first 11 dams, use data from previous time step to fill in missing data points
for i = 1:1:length(Dam)
    for k = 1:1:length(Plant.Data.Hydro.PowerGen(:,i))
        I = find(strcmp(Equipment{i+x},DamNames));
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

%Loading Temporary Electrical Demand
outputE = zeros(length(Dam),1);
for i = 1:1:length(Dam)
    Eff = Plant.Generator(i).VariableStruct.MaxGenCapacity/(Plant.Generator(i).VariableStruct.MaxGenFlow*Plant.Generator(i).VariableStruct.MaxHead/0.01181);%Power (kW)/ideal power in kW
    outputE(i) = Eff*Plant.Generator(i).VariableStruct.MaxHead*84.674;%Power (kW) = efficiency(%) * Head (ft) * 84.674 kJ/ (1000ft^3*ft)
end
Plant.Data.Demand.E = Plant.Data.Hydro.PowerFlow*outputE; %OutputPower = PowerFlow (1000 ft^3/s) * outputE
