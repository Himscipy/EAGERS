function WYHydroForecast(Date)
global Plant CurrentState DateSim TestData
hydroforecast = false;
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        hydroforecast = true;
    end
end
if isfield(Plant,'WYForecast') 
    nY = length(Plant.WYForecast);
    prevForcast = Plant.WYForecast{nY};
    if prevForcast.Timestamp(1)<=Date && prevForcast.Timestamp(end)>Date
        hydroforecast = false; %don't need to re-do
    end
end
if hydroforecast
    %first create the yearly dispatch data 
    % i.e. run Dispatch loop with updated information
    %these will be used for set points in the actual dispatch
    
    holdPlant = Plant; %save current plant so matrices do not have to be re-loaded
    holdTestData = TestData;%save to avoid re-interpolating at simulation frequency
    Plant.optimoptions.Horizon = 364*24; %Yearly Horizon
    Plant.optimoptions.Resolution = 7*24; %Week Resolution
    if isfield(Plant,'WYForecast')
        Plant = rmfield(Plant,'WYForecast');
    end
    D = datevec(Date);
    if D(2)<10
        year = D(1)-1;
    else
        year = D(1);
    end
    Date = datenum([year 10 1 1 0 0]);
    ForecastTime = Date+[0;buildTimeVector(Plant.optimoptions)/24];
    DateSim = Date(1);
    initializeOptimization
    if ~isfield(holdPlant,'WYForecast') || isempty(holdPlant.WYForecast)
        nY = 1;
        CurrentState.Hydro = zeros(1,length(holdPlant.subNet.Hydro.nodes));%SOC for reserviors,  IC is the initial power production
        for n=1:1:length(holdPlant.subNet.Hydro.nodes)
            i = holdPlant.subNet.Hydro.Equipment{n};
            if strcmp(holdPlant.Generator(i).Type,'Hydro Storage')
                CurrentState.Hydro(n) = holdPlant.Generator(i).VariableStruct.StartWYstate*holdPlant.Generator(i).QPform.Stor.UsableSize;
            end 
        end 
    else
        nY = length(holdPlant.WYForecast)+1;
    end
    reloadLast24hour(DateSim,Plant.optimoptions.Resolution)%re-load the previous 24 hours
    interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Horizon/24,0.00);%create test data at correct frequency
    Data = GetCurrentData(DateSim); 
    automaticInitialCondition(Data);
    Forecast = updateForecast(ForecastTime(2:end),true);
    if nY>1
        Solution = DispatchLoop(ForecastTime,Forecast,Plant.WYForecast{nY-1});
    else 
        Solution = DispatchLoop(ForecastTime,Forecast,[]);
    end 
    Solution.Timestamp = ForecastTime;
    Solution.hydroSOC = [CurrentState.Hydro;Solution.hydroSOC];
    disp(strcat('Water Year Forecast Completed for ',num2str(year),':',num2str(year+1)))
    Plant = holdPlant; %put original Plant with matrices back in place
    Plant.WYForecast{nY} = Solution;
    TestData = holdTestData;
end
end%Ends funtions WYHydroForecast