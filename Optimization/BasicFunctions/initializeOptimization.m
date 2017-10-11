global Plant
%% Load generators, & build QP matrices
%DateSim: Current time in the simulation.
%GenAvailTime: The next time a generator is available to turn on
%RestartTime: The amount of time necessary to wait for restarting after a generator has turned off.
if ~isfield(Plant.Data,'HistProf') || isempty(Plant.Data.HistProf)
    if isfield(Plant.Data,'Temperature')
        Plant.Data.HistProf.Temperature = TypicalDay([],Plant.Data.Timestamp,Plant.Data.Temperature);
    else
        Date = linspace(datenum([2017,1,1,1,0,0]),datenum([2018,1,1]),8760)'; %hours of 2017
        Plant.Data.HistProf.Temperature = TypicalDay([],Date,Plant.Weather.Tdb);
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
if strcmp(Plant.optimoptions.solver,'NREL')
    %do nothing
else
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
end