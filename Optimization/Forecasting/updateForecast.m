function Forecast = updateForecast(Date,Data)
%Date is the date number, Time is a vector of times (in hours)
global Plant
switch Plant.optimoptions.forecast
    case 'SNIWPE'
        Forecast = CreateSNIWPEForecast(Date,Data);
    case 'ARIMA'
        Forecast = CreateARIMAForecast(Date,Data);
    case 'NeuralNet'
        %%
    case 'Surface'
        Forecast = CreateForecast(Date,Data);
    case 'Perfect'
        Forecast = GetCurrentData(Date);
end

nS = length(Date);
nG = length(Plant.Generator);
Forecast.Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Forecast.Renewable(:,i) = RenewableOutput(Plant.Generator(i).VariableStruct,Date,'Predict');
    end
end