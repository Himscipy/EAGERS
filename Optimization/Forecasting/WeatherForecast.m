function Forecast = WeatherForecast(Date)
global Plant Last24hour
A = datevec(Last24hour.Timestamp(end));
hour = (Date-floor(Last24hour.Timestamp(end)))*24;
F = fieldnames(Last24hour.Weather);
for j = 1:1:length(F)
    YestFit = interp1(linspace(0,24,length(Last24hour.Weather.(F{j})))',Last24hour.Weather.(F{j}),mod((Date-Last24hour.Timestamp(end))*24,24));
    HistFit = interp1(0:24,[Plant.Data.HistProf.(F{j})(A(2),end),Plant.Data.HistProf.(F{j})(A(2),:)],mod(hour,24));
    W = interp1([Last24hour.Timestamp(end),Last24hour.Timestamp(end)+3,Last24hour.Timestamp(end)+100],[0.9,0,0],Date);%weight between yesterday and historical average
    W(isnan(W)) = 0;
    Forecast.(F{j}) =  (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
end
end %ends weatherForecast