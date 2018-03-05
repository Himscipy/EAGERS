function [Solution,Forecast] = SingleOptimization(Date,LastSolution)
% This function runs a single optimization over the horizon of the dispatch optimization
% Date is a vector of timestamps representing the points in time to be forecast and optimized
% LastSolution is optional. If it is not empty, [], it is a structure of the solution 1 step prior
% It returns Solution, a structure organizing the generator setpoints, line
% transfer flow, excess heat rejected, and hydro state of charge, and it
% returns the forecast that was optimized.
% If the Optimization has not been run previously it loads the optimization
% matrices and automatically determines an initial condition. If it has
% been run previously it uses the end point of a prior solution that is
% closest to the first date in ForecastTime as the initial condition. 
% This function then generates a forecast of the demands and optimizes
% according to that forecast.
global Plant DateSim CurrentState
DateSim = Date(1);
%It has not been run at this timestep (+- resolution)
if ~isfield(Plant,'Dispatch') || isempty(Plant.Dispatch) || ~isfield(Plant.Dispatch,'Timestamp') || ~any(Plant.Dispatch.Timestamp==(DateSim - Plant.optimoptions.Resolution/24))
    if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
        initializeOptimization
    end
    reloadLast24hour(DateSim,Plant.optimoptions.Resolution)%re-load the previous 24 hours
    interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Horizon/24,0.00);%create test data at correct frequency
    Data = updateForecast(DateSim);
    automaticInitialCondition(Data);
else
    t = max(linspace(1,length(Plant.Dispatch.Timestamp))'.*(Plant.Dispatch.Timestamp<DateSim)); %index preceeding current step
    CurrentState.Generators = Plant.Dispatch.GeneratorState(t,:);
end
Forecast = updateForecast(Date(2:end));
Solution = DispatchLoop(Date,Forecast,LastSolution);
end%ends Function single Optimization