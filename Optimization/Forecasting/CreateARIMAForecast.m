function Forecast = CreateARIMAForecast(Date,RealData)
global Plant Last24hour
Forecast.Timestamp = Date;
A = datevec(RealData.Timestamp);

%create forward forecast of temperature
hour = ([RealData.Timestamp; Date]-floor(RealData.Timestamp))*24;
YestFit = interp1(linspace(0,24,length(Last24hour.Temperature)+1)',[Last24hour.Temperature;RealData.Temperature],mod([0;(Date-RealData.Timestamp)*24],24));
HistFit = interp1(0:24,[Plant.Data.HistProf.Temperature(A(2),end),Plant.Data.HistProf.Temperature(A(2),:)],mod(hour,24));

W = interp1([RealData.Timestamp,RealData.Timestamp+3,RealData.Timestamp+100],[0.9,0,0],[RealData.Timestamp;Date]);%weight between yesterday and historical average
W(isnan(W)) = 0;
Tforecast = (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
S = fieldnames(Last24hour.Demand);
a = 60.9/100;
b= 45.3/100;
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [f,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.Demand.(S{s}) = zeros(length(Date),nD);
    r = Last24hour.Demand.(S{s});
    Date2 = [Last24hour.Timestamp;linspace(RealData.Timestamp,RealData.Timestamp+Plant.optimoptions.Horizon/24,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution+1)'];
    r(length(Date2),:) = 0;%pre-allocate
    d1 = r(2:end,:) - r(1:end-1,:);%kW/timestep
    for i = (length(Last24hour.Demand.(S{s})(:,1))+1):length(Date2)
        d1(i-1,:) = a*d1(i-2,:) + b*d1(i-f,:);%update the delta to include the new prediction 
        r(i,:) = d1(i-1,:) + r(i-1,:);
    end
    inter = interp1(Date2,r,Date);
    Forecast.Demand.(S{s}) = inter;
end
Forecast.Temperature = Tforecast(2:end);%exclude forecast for current time

if isfield(Plant.Data,'Hydro')
    fields = {'SpillFlow','OutFlow','InFlow','SourceSink'};
    for i = 1:1:length(fields)
        YestFit = interp1(linspace(0,24,length(Last24hour.Hydro.(fields{i}))+1)',[Last24hour.Hydro.(fields{i});RealData.Hydro.(fields{i})],mod([0;(Date-RealData.Timestamp)*24],24));
        HistFit = interp1(0:24,[Plant.Data.HistProf.Hydro.(fields{i})(A(2),end),Plant.Data.HistProf.Hydro.(fields{i})(A(2),:)],mod(hour,24));
        Forecast.Hydro.(fields{i}) = (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
    end
end