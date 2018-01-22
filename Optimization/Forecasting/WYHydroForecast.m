function WYHydroForecast(Date)
global Plant
hydroforecast = false;
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        hydroforecast = true;
    end
end
if hydroforecast
    %first create the yearly dispatch data 
    % i.e. run Dispatch loop with updated information
    %these will be used for set points in the actual dispatch
    if ~isfield(Plant,'WYForecast')
        Plant.WYForecast = []; %Start New Water Year
        nY = 1;
    else
        nY = length(Plant.WYForecast)+1;
    end
    OldOptions = Plant.optimoptions;
    Plant.optimoptions.Horizon = 364*24; %Yearly Horizon
    Plant.optimoptions.Resolution = 7*24; %Week Resolution
    Plant.subNet = []; %will force a re-construction of optimization matrices with current options
    D = datevec(Date);
    if D(2)<10
        Date = datenum([D(1)-1 10 1 1 0 0]);
    else
        Date = datenum([D(1) 10 1 1 0 0]);
    end
    ForecastTime = Date+[0;buildTimeVector(Plant.optimoptions)/24];
    if nY>1
        Solution = SingleOptimization(ForecastTime,Plant.WYForecast{nY-1});
    else 
        Solution = SingleOptimization(ForecastTime,[]);
    end 
    Solution.Dates = ForecastTime;
    Plant.WYForecast{nY} = Solution;
    Plant.optimoptions = OldOptions;
    Plant.subNet = []; %will force a re-construction of optimization matrices with current options
    initializeOptimization
    sprintf('WYForecast Completed')
end
end%Ends funtions WYHydroForecast