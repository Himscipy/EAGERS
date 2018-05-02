function [ForecastTime,Dispatch,HistoryTime,History] = SingleOptimizationNREL(Date)
global Plant TestData
n_g = length(Plant.Generator);%skip initialization
ic = zeros(1,n_g);
for i=1:1:n_g
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        ic(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
end
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    Buildings = Plant.Building;
else
    Buildings = [];
end
TestData.RealTimeData = interpolate_data(TestData,Plant.optimoptions.Resolution*3600,0.00);%create test data at correct frequency
ForecastTime = Date+[0;build_time_vector(Plant.optimoptions)/24];%linspace(Date,DateEnd)';would need to re-do optimization matrices for this time vector
[Forecast,Plant.Generator,~] = update_forecast(Plant.Generator,Buildings,[],[],Plant.optimoptions,ForecastTime(2:end));
scaleCost = update_cost(ForecastTime(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
Dispatch = solver_nrel(Plant.Generator,Plant.optimoptions,ic(1:n_g),22,Forecast,scaleCost);
%     dt = (ForecastTime(2:end) - ForecastTime(1:end-1))*24;
%     Cost = sum(Dispatch(2:end,1).*scaleCost(:,1).*dt) + sum(Dispatch(2:end,3).*scaleCost(:,3).*dt)/Plant.Generator(3).Output.Electricity(end);
History = [];
HistoryTime = [];
end %Ends function SingleOptimizationNREL