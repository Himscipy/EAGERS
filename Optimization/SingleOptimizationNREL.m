function [ForecastTime,Dispatch,HistoryTime,History] = SingleOptimizationNREL(Date)
global Plant Last24hour DateSim CurrentState
DateSim = Date;
nG = length(Plant.Generator);%skip initialization
IC = zeros(1,nG);
for i=1:1:nG
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        IC(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
end
CurrentState.Generators=IC(1:nG);
CurrentState.Building = 22.22;
Last24hour = [];%re-load the previous 24 hours
TimeYesterday = linspace(DateSim-1,DateSim,ceil(24/Plant.optimoptions.Resolution)+1)';
Last24hour = GetHistoricalData(TimeYesterday);
interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Horizon/24,0.00);%create test data at correct frequency
ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
Forecast = updateForecast(ForecastTime(2:end));
scaleCost = updateGeneratorCost(ForecastTime(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
[Dispatch,~] = NRELoptimization2(CurrentState.Generators,CurrentState.Building,Forecast,scaleCost);
%     dt = (ForecastTime(2:end) - ForecastTime(1:end-1))*24;
%     Cost = sum(Dispatch(2:end,1).*scaleCost(:,1).*dt) + sum(Dispatch(2:end,3).*scaleCost(:,3).*dt)/Plant.Generator(3).Output.Electricity(end);
History = [];
HistoryTime = [];
end %Ends function SingleOptimizationNREL