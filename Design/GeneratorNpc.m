function npc = GeneratorNpc(scale,Gen0,DesignDay,Years)
%NPCSIMULATE Calculate a Net Present Cost value for a given generator size.
global TestData Last24hour testSystems Plant DateSim GENINDEX SYSINDEX
Plant.Costs.Equipment(GENINDEX).Cost = testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost*scale;
Plant.Costs.Equipment(GENINDEX).OandM = testSystems(SYSINDEX).Costs.Equipment(GENINDEX).OandM*scale;
GenNew = updateComponentSpec(Gen0,'UB',Gen0.Size*scale);
F = fieldnames(GenNew);
for i = 1:1:length(F)
    Plant.Generator(GENINDEX).(F{i}) = GenNew.(F{i});
end
Plant.optimoptions.method = 'Planning';
Plant.optimoptions.forecast = 'Perfect';% Perfect forecast pulls directly from TestData
if isfield(Plant,'Building') &&  ~isempty(Plant.Building)
    Plant.optimoptions.forecast = 'Building';
end
Plant.optimoptions.Interval = floor(TestData.Timestamp(end)-TestData.Timestamp(1));
if DesignDay% If design days option is selected, optimize the 1st day of the month, and assume the rest of the month to be identical
    Plant.optimoptions.endSOC = 'Initial';% Constrain the final SOC of any storage device to be equal to theinitial charge so that days 2-30 of each month do not over-deplete storage.
    Plant.optimoptions.Horizon = max(24,Plant.optimoptions.Horizon);% make the horizon at least 1 day
    Plant.subNet = [];% empty the optimization matrices so that they are rebuilt with the new endSOC constraint
    DateSim = TestData.Timestamp(1);% set the starting date
    initializeOptimization% load optimizations
    NumSteps = preAllocateSpace('Design');% create Plant.Design structure with correct space
    interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Interval,0.00);% create test data at correct frequency
    STR = 'Optimizing Design Day Dispatch';
    DispatchWaitbar = waitbar(0, STR, 'Visible', 'on');
    Last24hour = [];% re-load the previous 24 hours
    TimeYesterday = linspace(DateSim-1, DateSim, ceil(24/Plant.optimoptions.Resolution)+1)';
    Last24hour = GetHistoricalData(TimeYesterday);
    automaticInitialCondition(GetCurrentData(DateSim));% set the initial conditions
    Plant.Design.Timestamp(1) = DateSim;
    ForecastTime = DateSim + [0; buildTimeVector(Plant.optimoptions)/24];
    Si = 1;
    while DateSim + Plant.optimoptions.Horizon/24 <= TestData.Timestamp(end)% loop to simulate days 1 to n in TestData
        D = datevec((DateSim));
        if Si == 1 || D(3) == 1% If it is the first step, or the first of the month run the actual optimization
            Forecast = updateForecast(ForecastTime(2:end));% function that creates demand vector with time intervals coresponding to those selected
            Solution = DispatchLoop(ForecastTime,Forecast,[]);
        else
            Forecast.Timestamp = ForecastTime(2:end);% otherwise just change the dates and use the preious solution
        end      
        StepDispatchForward(Si, ForecastTime, Forecast, Solution);% put solution into Plant.Design
        ForecastTime = round(864000*(ForecastTime+Plant.optimoptions.Horizon/24)) / 864000;% count forward by length of the horizon, rounded to nearest second
        Si = Si + length(ForecastTime)-1;
        waitbar(Si/NumSteps, DispatchWaitbar,strcat('Running Design Day Dispatch'));
    end
    close(DispatchWaitbar)
else
    Plant.optimoptions.endSOC = 'Flexible';% remove constraint on final SOC of storage
    Plant.subNet = [];%empty the optimization matrices so that they are rebuilt without endSOC constraint
    initializeOptimization
    preAllocateSpace('Design')
    if ~isfield(Plant,'Design') || isempty(Plant.Design)|| any(Plant.Design.Timestamp==0)
        % at least some points have not been run
        STR = 'Optimizing Dispatch Throughout Entire Year';
        DispatchWaitbar = waitbar(0, STR, 'Visible', 'on');
        RunSimulation(TestData.Timestamp(1),[]);
        close(DispatchWaitbar)
    end
end
Plant.optimoptions = testSystems(SYSINDEX).optimoptions;
testSystems(SYSINDEX).Design = Plant.Design;
F = fieldnames(GenNew);
for i = 1:1:length(F)
    testSystems(SYSINDEX).Generator(GENINDEX).(F{i}) = GenNew.(F{i});
end
DesignCosts(SYSINDEX,Years,Plant.Costs.Equipment);% update the costs, monthly costs & NPC for system k
npc = testSystems(SYSINDEX).Costs.NPC;
end%Ends function GeneratorNpc