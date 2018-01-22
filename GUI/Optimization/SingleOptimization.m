function Solution = SingleOptimization(ForecastTime,LastDispatch)
global Plant DateSim CurrentState RealTimeData
DateSim = ForecastTime(1);
%It has not been run at this timestep (+- resolution)
if ~isfield(Plant,'Dispatch') || isempty(Plant.Dispatch) || ~isfield(Plant.Dispatch,'Timestamp') || ~any(Plant.Dispatch.Timestamp==(DateSim - Plant.optimoptions.Resolution/24))
    if isempty(Plant.subNet)
        initializeOptimization
    end
    loadTestData
    resetLast24Hour
    timestamp = (DateSim:Plant.optimoptions.Resolution/24:(DateSim+Plant.optimoptions.Horizon/24))';
    RealTimeData = interpolateData(timestamp,0.00);%create test data at correct frequency
    if isempty(Plant.Building)
        Data = GetHistoricalData(DateSim); 
    else
        Data = updateForecast(DateSim,[]);
    end
    automaticInitialCondition(Data);
else
    t = max(linspace(1,length(Plant.Dispatch.Timestamp))'.*(Plant.Dispatch.Timestamp<DateSim)); %index preceeding current step
    CurrentState.Generators = Plant.Dispatch.GeneratorState(t,:);
    CurrentState.Lines = Plant.Dispatch.LineFlows(t,:); 
    CurrentState.Buildings = Plant.Dispatch.Buildings(t,:);
end
Forecast = updateForecast(ForecastTime(2:end),Data);
Solution = DispatchLoop(ForecastTime,Forecast,LastDispatch);
HeatRecovery = CalcHeatRecovery(Solution.Dispatch);
end%ends Function single Optimization