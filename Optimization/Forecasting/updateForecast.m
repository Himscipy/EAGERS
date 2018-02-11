function Forecast = updateForecast(Date)
%Date is the date number, Time is a vector of times (in hours)
global Plant
switch Plant.optimoptions.forecast
    case 'SNIWPE'
        Forecast = SNIWPEForecast(Date);
    case 'ARIMA'
        Forecast = ARIMAForecast(Date);
    case 'NeuralNet'
        %
    case 'Surface'
        Weather = WeatherForecast(Date);
        Forecast = SurfaceForecast(Date,Weather.Tdb);
    case 'Perfect'
        Forecast = GetCurrentData(Date);
    case 'Building'
        Weather = WeatherForecast(Date);
        nB = length(Plant.Building);
        for i = 1:1:nB
            Build = Plant.Building(i);
            Location = Plant.subNet.Electrical.Location(Build.QPform.Electrical.subnetNode);
            [Forecast.Building.InternalGains(:,i),Forecast.Building.ExternalGains(:,i),Equipment,InteriorLighting, ExteriorLighting, OtherLoads] = BuildingLoads(Build,Weather.irradDireNorm,Weather.irradDiffHorz,Location,Date);
            Forecast.Building.NonHVACelectric(:,i) = Equipment + InteriorLighting + ExteriorLighting + OtherLoads;
        end
end
Forecast.Timestamp = Date;
if isfield(Date,'Weather')
    Forecast.Weather = WeatherForecast(Date);
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
end
if isfield(Forecast,'Building')
    Forecast.Building.ExternalGains = ForecastExternalGains(Forecast);
    [Forecast.Building,Forecast.SRtarget] = ForecastBuilding(Weather,Date,Forecast.Building);
end
if isfield(Forecast,'Weather') && isfield(Forecast.Weather,'irradDireNorm')
    Forecast.Renewable = RenewableOutput(Date,Forecast.Weather.irradDireNorm);
end
if isfield(Plant,'subNet') && isfield(Plant.subNet,'Hydro')
    Forecast.Hydro.InFlow = ForecastHydro(Date,Forecast.Hydro.SourceSink);
end

if Plant.optimoptions.SpinReserve
    if isfield(Forecast,'Demand')
        Forecast.SRtarget = Forecast.SRtarget + Plant.optimoptions.SpinReservePerc/100*sum(Forecast.Demand.E,2);
    end
end
end%Ends function updateForecast