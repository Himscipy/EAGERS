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
Last24hour = [];
TimeYesterday = linspace(DateSim-1+Plant.optimoptions.Resolution/24,DateSim,24/Plant.optimoptions.Resolution)';
if Plant.Data.Timestamp(1)<(DateSim-1) && Plant.Data.Timestamp(end)>DateSim
    Last24hour = GetHistoricalData(TimeYesterday);
else %need to have this in terms of the first timestep
    Last24hour = GetHistoricalData(TimeYesterday+1);
    Last24hour.Timestamp = TimeYesterday;
end
Data = GetHistoricalData(DateSim); 
ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
Forecast = updateForecast(ForecastTime(2:end),Data);
scaleCost = updateGeneratorCost(ForecastTime(2:end)); %% All feedstock costs were assumed to be 1 when building matrices 
[Dispatch,~] = NRELoptimization2(CurrentState.Generators,CurrentState.Building,Forecast,scaleCost);
%     dt = (ForecastTime(2:end) - ForecastTime(1:end-1))*24;
%     Cost = sum(Dispatch(2:end,1).*scaleCost(:,1).*dt) + sum(Dispatch(2:end,3).*scaleCost(:,3).*dt)/Plant.Generator(3).Output.Electricity(end);
History = [];
HistoryTime = [];
end %Ends function SingleOptimizationNREL