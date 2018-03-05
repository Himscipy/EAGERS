function Forecast = updateForecast(varargin)
%Date is the date number, Time is a vector of times (in hours)
global Plant
Date = varargin{1};
if length(varargin)==1 || varargin{2} == false
    WYHydroForecast(Date(end));%if october 1st,Run a yearly forecast for hydrology
end
switch Plant.optimoptions.forecast
    case 'SNIWPE'
        Forecast = SNIWPEForecast(Date);
        Forecast.Weather = WeatherForecast(Date);
    case 'ARIMA'
        Forecast = ARIMAForecast(Date);
        Forecast.Weather = WeatherForecast(Date);
    case 'NeuralNet'
        %
    case 'Surface'
        Weather = WeatherForecast(Date);
        Forecast = SurfaceForecast(Date,Weather.Tdb);
        Forecast.Weather = Weather;
    case 'Perfect'
        Forecast = GetCurrentData(Date);
    case 'Building'
        Forecast.Timestamp = Date;
        Forecast.Weather = WeatherForecast(Date);
        nB = length(Plant.Building);
        for i = 1:1:nB
            Build = Plant.Building(i);
            Location = Plant.subNet.Electrical.Location(Build.QPform.Electrical.subnetNode);
            SG = SolarGain(Build,Date,Location,Forecast.Weather);
            Forecast.Building.ExternalGains = SG.Walls + SG.Roof;
            B_loads = BuildingLoads(Build,Date,SG);
            Forecast.Building.InternalGains(:,i) = B_loads.InternalGains;
            Forecast.Building.NonHVACelectric(:,i) = B_loads.Equipment + B_loads.InteriorLighting + B_loads.ExteriorLighting + B_loads.OtherLoads;
            if isfield(B_loads,'DCloads')
                Forecast.Building.DCloads = B_loads.DCloads;
            end
        end
end
if isfield(Forecast,'Building') && ~isfield(Forecast.Building,'ExternalGains')
    %needed when forecast method is something other than 'Building', and TestData does not have stored ExternalGains
    nB = length(Plant.Building);
    Forecast.Building.ExternalGains = zeros(length(Forecast.Timestamp),nB);
    for i = 1:1:nB
        Build = Plant.Building(i);
        Location = Plant.Building(i).QPform.Location;
        SG = SolarGain(Build,Forecast.Timestamp,Location,Forecast.Weather);
        Forecast.Building.ExternalGains(:,i) = SG.Walls + SG.Roof;
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
end
if isfield(Forecast,'Building') && ~strcmp(Plant.optimoptions.solver,'NREL')
    %needed for WSU optimization approach only
    [Forecast.Building,Forecast.SRtarget] = ForecastBuilding(Forecast.Weather,Date,Forecast.Building);
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