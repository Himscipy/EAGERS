global Plant Model_dir
Plant = [];
Plant.Name = 'ColumbiaBasin';
Plant.Data.Timestamp(:,1)= linspace(datenum([2015 10 2 1 0 0]),datenum([2016 10 2 0 0 0]),8784).';
Plant.Data.Temperature = 20*ones(length(Plant.Data.Timestamp),1);
Plant.Data.Holidays = [];
Plant.Data.HistProf = [];

%% Network Setup
%Brownlee, Oxbow, and Hells Canyon Dams ommitted

NodeNames = {'Plant1';'Plant2';'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa'};
         
DownRiver = {{'Plant2'};{''};'Bridgeport, Wa';'Azwell, Wa';''};
           
UpRiver =  {{''};{'Plant1'}};

Equipment = {'Plant1';'Plant2';'Grand Coulee';'Chief Joseph';'Wells'};

HydroConnection = {'Plant1';'Plant2';'Plant2'};

InstreamFlow = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10;];
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
Dam(3).VariableStruct.MaxGenCapacity = 851400; %power in kW
Dam(3).VariableStruct.RampUp =599318; %power change in kW/hr
Dam(3).VariableStruct.RampDown =613318; 
Dam(3).VariableStruct.MaxGenFlow = 219.2; %flow in 1000 cfs
Dam(3).VariableStruct.MaxSpillFlow = 185.5; %flow in 1000 cfs
Dam(3).VariableStruct.MaxHead = 76.65; 
Dam(3).VariableStruct.MinHead = 46.91; %guess minimum height allowable

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
load(fullfile(Model_dir,'Data','Hydro','allHistHydroGenData_2007_2016')); %Historical Generation Data(hourly)(kW)
% load(fullfile(Model_dir,'Data','Hydro','SourcesandSinks_2007_2016')); %Source and Sink Data; UpstreamOutFlow(t-T)-DownstreamInFlow(t)
%%columns 1->18 (kcfs or KW): Grand Coulee, Chief Joseph, Wells, Rocky Ridge, Rock Island, Wanapum, PriestRiver, McNary, John Day, The Dalles, Bonneville, Brownlee, Oxbow, Hells Canyon, Lower Granite, Little Goose, Lower Monumental, Ice Harbor;
%SourceSink = Mass balance; Sinks and Sources in river segments; Negative values = sinks; Positive values = sources;

%These are corrected data 
load(fullfile(Model_dir,'Data','Hydro','SourceSinkHrly'))
load(fullfile(Model_dir,'Data','Hydro','InFlowCrctd'))
load(fullfile(Model_dir,'Data','Hydro','OutFlowNMD'))
load(fullfile(Model_dir,'Data','Hydro','powerFlowNMD'))

for i = 1:1:length(Dam) 
    if exist(fullfile(Model_dir,'Data','Hydro',strcat(DamNames{i},'_','2007','_','2016','.mat')),'file')
        name = fullfile(Model_dir,'Data','Hydro',strcat(DamNames{i},'_','2007','_','2016'));
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
Plant.Network(2).Electrical.Load = 1; %put all the load at the first node

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
