function [forecast,gen,building] = update_forecast(gen,building,cool_tower,subnet,options,date,hist_prof,prev_data,now_data,future_data)
%Date is the date number, Time is a vector of times (in hours)
[forecast.WYForecast,gen] = water_year_forecast(gen,building,cool_tower,subnet,options,date(end),hist_prof,prev_data,now_data,future_data);%if october 1st,Run a yearly forecast for hydrology
switch options.forecast
    case 'ARMA'
        forecast = arma(date,prev_data);
        forecast.Weather = weather_forecast(prev_data,hist_prof,date);
    case 'ARIMA'
        forecast = arima_wsu(date,prev_data,options);
        forecast.Weather = weather_forecast(prev_data,hist_prof,date);
    case 'NeuralNet'
        %
    case 'Surface'
        Weather = weather_forecast(prev_data,hist_prof,date);
        forecast = surface_forecast(prev_data,hist_prof.Demand,date,Weather.Tdb,[]);
        forecast.Weather = Weather;
    case 'Perfect'
        if length(date) == 1
            forecast = now_data;
        else
            forecast = future_data;
        end
    case 'Building'
        forecast.Timestamp = date;
        forecast.Weather = weather_forecast(prev_data,hist_prof,date);
        n_b = length(building);
        for i = 1:1:n_b
            sgain = solar_gain(building(i),date,building(i).QPform.Location,forecast.Weather);
            forecast.Building.ExternalGains = sgain.Walls + sgain.Roof;
            b_loads = building_loads(building(i),date,sgain);
            forecast.Building.InternalGains(:,i) = b_loads.InternalGains;
            forecast.Building.NonHVACelectric(:,i) = b_loads.Equipment + b_loads.InteriorLighting + b_loads.ExteriorLighting + b_loads.OtherLoads;
            if isfield(b_loads,'DCloads')
                forecast.Building.DCloads = b_loads.DCloads;
            end
        end
end
if isfield(forecast,'Building') && ~isfield(forecast.Building,'ExternalGains')
    %needed when forecast method is something other than 'Building', and TestData does not have stored ExternalGains
    n_b = length(building);
    forecast.Building.ExternalGains = zeros(length(forecast.Timestamp),n_b);
    for i = 1:1:n_b
        sgain = solar_gain(building(i),forecast.Timestamp,building(i).QPform.Location,forecast.Weather);
        forecast.Building.ExternalGains(:,i) = sgain.Walls + sgain.Roof;
    end
end
if strcmp(options.method,'Dispatch')
    f = fieldnames(now_data);
    for j = 1:1:length(f)
        if isstruct(now_data.(f{j}))
            S = fieldnames(now_data.(f{j}));
            for i = 1:1:length(S)
                forecast.(f{j}).(S{i})(1,:) = now_data.(f{j}).(S{i});
            end
        else
            forecast.(f{j})(1,:) = now_data.(f{j});
        end
    end
end
if isfield(forecast,'Building') && ~strcmp(options.solver,'NREL')
    %needed for WSU optimization approach only
    %%Need warm-up period if not currently running the model
    for i = 1:1:n_b
        if  abs(round(864000*(building(i).Timestamp+options.Resolution/24))/864000 - date(1))>1e-5
            if length(date) == 1
                wu_date = linspace(date(1) + options.Resolution/24,date(1)+1,24)';
            else
                wu_date = linspace(date(1),date(1)+1 - options.Resolution/24,24)';
            end
            wu_weather = weather_forecast(prev_data,hist_prof,wu_date);
            building(i) = building_warmup(building(i),wu_weather,options.Resolution/24,wu_date,6);
        end
    end
    forecast.Building = forecast_building(forecast.Weather,date,forecast.Building,building);
end
if isfield(forecast,'Weather') && isfield(forecast.Weather,'irradDireNorm')
    forecast.Renewable = renewable_output(gen,subnet,date,forecast.Weather.irradDireNorm);
end
if isfield(subnet,'Hydro')
    forecast.Hydro.InFlow = forecast_hydro(prev_data,date,forecast.Hydro.SourceSink,subnet,options.Resolution/24);
end

if options.SpinReserve
    if isfield(forecast,'Demand')
        forecast.SRtarget = forecast.SRtarget + options.SpinReservePerc/100*sum(forecast.Demand.E,2);
    end
end
end%Ends function update_forecast