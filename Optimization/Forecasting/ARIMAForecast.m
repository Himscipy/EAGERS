function Forecast = ARIMAForecast(Date)
global Plant Last24hour
S = fieldnames(Last24hour.Demand);
a = 60.9/100;
b= 45.3/100;
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [f,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.Demand.(S{s}) = zeros(length(Date),nD);
    r = Last24hour.Demand.(S{s});
    Date2 = [Last24hour.Timestamp;linspace(Last24hour.Timestamp(end)+Plant.optimoptions.Resolution/24,Last24hour.Timestamp(end)+Plant.optimoptions.Horizon/24,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution+1)'];
    r(length(Date2),:) = 0;%pre-allocate
    d1 = r(2:end,:) - r(1:end-1,:);%kW/timestep
    for i = (length(Last24hour.Demand.(S{s})(:,1))+1):length(Date2)
        d1(i-1,:) = a*d1(i-2,:) + b*d1(i-f,:);%update the delta to include the new prediction 
        r(i,:) = d1(i-1,:) + r(i-1,:);
    end
    inter = interp1(Date2,r,Date);
    Forecast.Demand.(S{s}) = inter;
end