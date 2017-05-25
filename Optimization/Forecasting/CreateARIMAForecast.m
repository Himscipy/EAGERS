function Forecast = CreateARIMAForecast(Date,RealData)
global Plant Last24hour
Forecast.Timestamp = Date;
Days = ceil(Date(end)-floor(RealData.Timestamp)); %this says that you may need to use surface fits for multiple days.
s12hour = nnz(Date<(RealData.Timestamp+.5)); %steps in next 12 hours
A = datevec(RealData.Timestamp);
dateDays = floor(RealData.Timestamp)-1:ceil(RealData.Timestamp+Days);
Weekday = zeros(length(dateDays),1);
for i = 1:1:length(dateDays)
    % Is dateDays(i) a weekday and not a holiday?
    if weekday(dateDays(i))<=6 && weekday(dateDays(i))>=2 && ~any(Plant.Data.Holidays==dateDays(i))
        Weekday(i) = 1;
    else
        Weekday(i) = 0;
    end
    % Is the historical data from a year in the past? (DOES NOT HANDLE
    % HISTORICAL DATA FROM MORE THAN 1 YEAR IN THE PAST.)
    if Weekday(i)==1 && any(Plant.Data.Holidays == dateDays(i)+365) && dateDays(i)<min(Plant.Data.Holidays)
        Weekday(i) = 0;
    end
end
monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};

%create forward forecast of temperature
hour = ([RealData.Timestamp; Date]-floor(RealData.Timestamp))*24;
YestFit = interp1(linspace(0,24,length(Last24hour.Temperature)+1)',[Last24hour.Temperature;RealData.Temperature],mod([0;(Date-RealData.Timestamp)*24],24));
HistFit = interp1(0:24,[Plant.Data.HistProf.Temperature(A(2),end),Plant.Data.HistProf.Temperature(A(2),:)],mod(hour,24));

W = interp1([RealData.Timestamp,RealData.Timestamp+3,RealData.Timestamp+100],[0.9,0,0],[RealData.Timestamp;Date]);%weight between yesterday and historical average
W(isnan(W)) = 0;
Tforecast = (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
S = fieldnames(Last24hour.Demand);
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [f,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.Demand.(S{s}) = zeros(length(Date),nD);
    for k = 1:1:nD
        D = Last24hour.Demand.(S{s});
        d1 = diff(D(1:end));%kW/timestep
        r=[];
        a = 60.9;
        b= 45.3;
        h0 = Date(end)-Date(1);
        h = round(h0*24/Plant.optimoptions.Resolution);
        start = f+1;
        finish = f + h;
        for i = start:finish%For next 24 hours
            pred = a*.01*d1(i-2) + b*.01*d1(i-f);%1 timestep and 24 hours
            if isempty(r)
                r(:,end+1) = pred + D(i-1);%If 1st prediction, add predicted delta to real E @ hr 24
            else
                r(:,end+1) = pred + r(end);%If not 1st, add predicted delta to predicted E in last timestep
            end
            d1(end+1,:) = pred; %update the delta to include the new prediction 
        end
        Date2 = linspace(Last24hour.Timestamp(1) +1,Last24hour.Timestamp(end)+h0,h);
        inter = interp1(Date2,r,Date);
        Forecast.Demand.(S{s}) = inter';
    end
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