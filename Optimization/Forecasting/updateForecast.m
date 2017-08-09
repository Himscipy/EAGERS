function Forecast = updateForecast(Date,Data)
%Date is the date number, Time is a vector of times (in hours)
global Plant Last24hour
nS = length(Date);
nG = length(Plant.Generator);
nB = length(Plant.Building);
if ismember(Plant.optimoptions.forecast,{'Surface';'ARIMA';'Building'})%create forward forecast of temperature
    if ~isempty(Data) && isfield(Last24hour,'Temperature')
        A = datevec(Data.Timestamp);
        hour = ([Data.Timestamp; Date]-floor(Data.Timestamp))*24;
        YestFit = interp1(linspace(0,24,length(Last24hour.Temperature)+1)',[Last24hour.Temperature;Data.Temperature],mod([0;(Date-Data.Timestamp)*24],24));
        HistFit = interp1(0:24,[Plant.Data.HistProf.Temperature(A(2),end),Plant.Data.HistProf.Temperature(A(2),:)],mod(hour,24));
        W = interp1([Data.Timestamp,Data.Timestamp+3,Data.Timestamp+100],[0.9,0,0],[Data.Timestamp;Date]);%weight between yesterday and historical average
        W(isnan(W)) = 0;
        Tforecast =  (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
    else
        D = datevec(Date(1));
        h_of_y = 24*(Date - datenum([D(1),1,1])); %hour since start of year
        Tforecast = interp1((0:1:8760)',[Plant.Weather.Tdb(1);Plant.Weather.Tdb],h_of_y);
        if length(h_of_y) ==1
            Tforecast(2,1) = Tforecast;
        end
    end
end
switch Plant.optimoptions.forecast
    case 'SNIWPE'
        Forecast = SNIWPEForecast(Date,Data);
    case 'ARIMA'
        Forecast = ARIMAForecast(Date,Data);
    case 'NeuralNet'
        %%
    case 'Surface'
        Forecast = SurfaceForecast(Date,Data,Tforecast);
    case 'Perfect'
        Forecast = GetCurrentData(Date);
    case 'Building'
        for i = 1:1:nB
            Forecast.Building(i) = ForecastBuilding(Plant.Building(i),Plant.Weather,Date,Plant.Building(i).QPform.Location);
        end
end
Forecast.Timestamp = Date;
if ismember(Plant.optimoptions.forecast,{'Surface';'ARIMA';'Building'})%create forward forecast of temperature
    Forecast.Temperature = Tforecast(2:end);%exclude forecast for current time
end
Forecast.Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Forecast.Renewable(:,i) = RenewableOutput(i,Date,'Predict');
    end
end
if isfield(Plant.subNet,'Hydro')
    Forecast.Hydro.InFlow = ForecastHydro(Date,Data.Timestamp);
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