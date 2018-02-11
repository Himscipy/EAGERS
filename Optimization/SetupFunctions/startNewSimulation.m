global Plant  DateSim Last24hour CurrentState
startConstraints
if ~isempty(Plant.Dispatch.hydroSOC(1,:))
    WYHydroForecast(DateSim); %Water Year Forecast for Hydro
else
    if any(strcmp(Plant.optimoptions.method,{'Dispatch';'Planning'}))
        interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Interval,0.00);%create test data at correct frequency
    else %actual control
        interpolateData(Plant.optimoptions.Tmpc,Plant.optimoptions.Interval,0.00);%create test data at correct frequency
    end
end
if Compare%%If running for comparisson of mixed and non-mixed integer
    Plant.Predicted.GenDispcQP = zeros(size(Plant.Predicted.GenDisp));
    Plant.Predicted.CostcQP = zeros(NumSteps,1);
    timers_cQP = zeros(NumSteps,3);
end
Data = GetCurrentData(DateSim);
Last24hour = [];%re-load the previous 24 hours
TimeYesterday = linspace(DateSim-1,DateSim,ceil(24/Plant.optimoptions.Resolution)+1)';
Last24hour = GetHistoricalData(TimeYesterday);
% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;
if K ==1
    manualInitialCondition;
else
    automaticInitialCondition(Data);
end
Plant.Dispatch.GeneratorState(1,:) = CurrentState.Generators;
Plant.Dispatch.Buildings(1,:) = CurrentState.Buildings(1,:);
Plant.Dispatch.Timestamp(1) = DateSim;