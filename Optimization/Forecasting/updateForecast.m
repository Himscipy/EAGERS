function Forecast = updateForecast(Date,Data)
%Date is the date number, Time is a vector of times (in hours)
global Plant
nS = length(Date);
nG = length(Plant.Generator);
nB = length(Plant.Building);
switch Plant.optimoptions.forecast
    case 'SNIWPE'
        Forecast = SNIWPEForecast(Date,Data);
    case 'ARIMA'
        Forecast = ARIMAForecast(Date,Data);
    case 'NeuralNet'
        %
    case 'Surface'
        Forecast = SurfaceForecast(Date,WeatherForecast(Date,'Tdb'));
    case 'Perfect'
        Forecast = GetCurrentData(Date);
    case 'Building'
        for i = 1:1:nB
            Forecast.Building(i) = ForecastBuilding(Plant.Building(i),Plant.Weather,Date,Plant.Building(i).QPform.Location);
        end
end
Forecast.Timestamp = Date;
% Forecast.Temperature = WeatherForecast(Date,'Temperature');
S = {'Tdb';'Twb';'irradDireNorm'};
for j = 1:1:length(S)
    Forecast.Weather.(S{j}) = WeatherForecast(Date,S{j});
end
Forecast.Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Forecast.Renewable(:,i) = RenewableOutput(i,Date,'Predict',Forecast.Weather.irradDireNorm);
    end
end
if strcmp(Plant.optimoptions.method,'Dispatch')
    Data = GetCurrentData(Date(1));
    F = fieldnames(Data);
    for j = 1:1:length(F)
        if isstruct(Data.(F{j}))
            S = fieldnames(Data.(F{j}));
            for i = 1:1:length(S)
                Forecast.(F{j}).(S{i})(1,:) = Data.(F{j}).(S{i});
            end
        else
            Forecast.(F{j})(1,:) = Data.(F{j});
        end
    end
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Forecast.Renewable(1,i) = RenewableOutput(i,Date(1),'Actual',Forecast.Weather.irradDireNorm(1));
        end
    end
end

if isfield(Plant.subNet,'Hydro')
    Forecast.Hydro.InFlow = ForecastHydro(Date,Forecast.Hydro.SourceSink);
end

if Plant.optimoptions.SpinReserve
    if isfield(Forecast,'Demand')
        Forecast.SRtarget = Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E,2);
    else
        Forecast.SRtarget  = Forecast.Building(1).E0;
        for i = 2:1:nB
            Forecast.SRtarget = Forecast.SRtarget + Forecast.Building(i).E0;
        end
    end
end
end%Ends function updateForecast