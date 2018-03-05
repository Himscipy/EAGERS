global Plant  DateSim CurrentState
startConstraints
if any(strcmp(Plant.optimoptions.method,{'Dispatch';'Planning'}))
    interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Interval,0.00);%create test data at correct frequency
else %actual control
    interpolateData(Plant.optimoptions.Tmpc,Plant.optimoptions.Interval,0.00);%create test data at correct frequency
end
Data = updateForecast(DateSim);%Data = GetCurrentData(DateSim);
reloadLast24hour(DateSim,Plant.optimoptions.Resolution)%re-load the previous 24 hours
% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;
if K ==1
    manualInitialCondition;
else
    
    if strcmp(Plant.optimoptions.solver,'ANN')
        Plant.optimoptions.solver = 'quadprog';
        automaticInitialCondition(Data);
        Plant.optimoptions.solver = 'ANN';
    else
        automaticInitialCondition(Data);
    end
end
Plant.Dispatch.GeneratorState(1,:) = CurrentState.Generators;
if isfield(CurrentState,'Buildings')
    Plant.Dispatch.Buildings(1,:) = CurrentState.Buildings(1,:);
end
if isfield(CurrentState,'Hydro')
    Plant.Dispatch.hydroSOC(1,:) = CurrentState.Hydro(1,:);
end
Plant.Dispatch.Timestamp(1) = DateSim;