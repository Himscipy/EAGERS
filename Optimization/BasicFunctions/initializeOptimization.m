global Plant
%% Load generators, & build QP matrices
%DateSim: Current time in the simulation.
%GenAvailTime: The next time a generator is available to turn on
%RestartTime: The amount of time necessary to wait for restarting after a generator has turned off.
if strcmp(Plant.optimoptions.forecast,'Surface') && (~isfield(Plant.Data,'HistProf') || isempty(Plant.Data.HistProf))
    calculateHistoricalFit %% calculate surface fits used in forecasting
end

loadGenerator % Loads generators & build optimization matrices, stored in Plant.Generator
build_subNet
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
Plant.OpMatA = buildMatrices('OpMatA',dt); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('OpMatB',dt);%build quadratic programming matrices for FitB
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