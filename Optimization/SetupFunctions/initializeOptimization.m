global Plant Model_dir
allFieldNames = {'Name';'Data';'Generator';'Building';'Weather';'Network';'Costs';'optimoptions';'subNet';'OpMatA';'OpMatB';'OneStep';'Online';'Design';'Dispatch';'Predicted';'RunData';'Baseline';'Market'};
fNames = fieldnames(Plant);
for i = 1:1:length(allFieldNames)
    if ~any(strcmp(allFieldNames{i},fNames))
        Plant.(allFieldNames{i}) = [];
    end
end
if isempty(Plant.Weather)
    load(fullfile(Model_dir,'System Library','Weather','4A.mat'))
    Plant.Weather = weather;
end

if ~isfield(Plant.Data,'Weather') || ~isfield(Plant.Data.Weather,'Tdb') %% if there was not weather data lined up with the historical data, add it
    D = datevec(Plant.Data.Timestamp(1));
    h_of_y = linspace(0,8760,length(Plant.Weather.Tdb)+1)';% Hour since start of year
    t = mod(24*(Plant.Data.Timestamp - datenum([D(1),1,1])),8760);% Hour since start of year
    Plant.Data.Weather.Tdb = interp1(h_of_y,[Plant.Weather.Tdb(1); Plant.Weather.Tdb],t); 
end

if ~isfield(Plant.Data,'HistProf') || isempty(Plant.Data.HistProf)
    S = {'Tdb';'Twb';'irradDireNorm';};
    for j = 1:1:length(S)
        if isfield(Plant.Data,'Weather') && isfield(Plant.Data.Weather,S{j})
            Plant.Data.HistProf.(S{j}) = TypicalDay([],Plant.Data.Timestamp,Plant.Data.Weather.(S{j}));
        else
            Date = linspace(datenum([2017,1,1,1,0,0]),datenum([2018,1,1]),8760)'; %hours of 2017
            Plant.Data.HistProf.(S{j}) = TypicalDay([],Date,Plant.Weather.(S{j}));
        end
    end
    if strcmp(Plant.optimoptions.forecast,'Surface') && isfield(Plant.Data,'Demand') && ~isempty(Plant.Data.Demand)
        F = fieldnames(Plant.Data.Demand);
        for i = 1:1:length(F)
            calculateHistoricalFit(F{i}); %% calculate surface fits used in forecasting
        end
    end
end
if isfield(Plant.Data,'Hydro')
    %spill, outflow, inflow, and source/sink by node
    nodes = length(Plant.Data.Hydro.SourceSink(1,:));
    for n = 1:1:nodes
%         Plant.Data.HistProf.Hydro.SpillFlow(n) = {TypicalDay([],Plant.Data.Hydro.Timestamp,Plant.Data.Hydro.SpillFlow)};
%         Plant.Data.HistProf.Hydro.OutFlow(n) = {TypicalDay([],Plant.Data.Hydro.Timestamp,Plant.Data.Hydro.OutFlow)};
%         Plant.Data.HistProf.Hydro.InFlow(n) = {TypicalDay([],Plant.Data.Hydro.Timestamp,Plant.Data.Hydro.InFlow)};
        Plant.Data.HistProf.Hydro.SourceSink(n) = {TypicalDay([],Plant.Data.Hydro.Timestamp,Plant.Data.Hydro.SourceSink(:,n))};
    end
end
%% Load generators, & build QP matrices
loadGenerator % Loads generators stored in Plant.Generator
loadBuilding %Loads any buildings into QPform
build_subNet
findbuffer %need to wait until everything else has been loaded to load the buffer, because it is reliant on how quickly everything else can charge the system
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
Plant.OpMatA = buildMatrices('A',dt); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B',dt);%build quadratic programming matrices for FitB
Plant.OneStep = buildMatrices1Step;%build quadratic programming matrices for 1Step at interval spacing of dt

if strcmp(Plant.optimoptions.method,'Control')
    A.Horizon = Plant.optimoptions.Resolution;%the horizon is the resolution
    A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = buildTimeVector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        Plant.Online(t) = buildMatrices('OpMatB',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end

