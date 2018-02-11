function TestWeather = interpolateWeather(weather,Date)
TestWeather.Name = weather.Name;
if ~isfield(weather,'Timestamp')%Assume it is a typical meteorological year
    weather.Timestamp = linspace(datenum([2017,1,1,8760/length(weather.Tdb),0,0]),datenum([2018,1,1]),length(weather.Tdb))';
end
nS = length(Date);
S = fieldnames(weather);
S = S(~strcmp('Timestamp',S));
for j = 1:1:length(S)
    if isnumeric(weather.(S{j}))
        TestWeather.(S{j}) = nan(nS,1);
    end
end
if Date(1)>=weather.Timestamp(1) && Date(end)<=weather.Timestamp(end)
    for j = 1:1:length(S)
        if isnumeric(weather.(S{j}))
            TestWeather.(S{j}) = interp1(weather.Timestamp,weather.(S{j}),Date)'; 
        end
    end
else %need to move solar data timestamp to line up with TestData timestamp
    D1 = datevec(Date(1));
    D2 = datevec(weather.Timestamp(1));
    
    weather.Timestamp = weather.Timestamp + datenum([D1(1),1,1]) - datenum([D2(1),1,1]);%days between start of years for testdata and solar data respectively
    if datenum([D1(1),D2(2),D2(3),D2(4),D2(5),D2(6)])>Date(1)
        weather.Timestamp = weather.Timestamp - (datenum([D1(1),1,1]) - datenum([D1(1)-1,1,1]));
    end
    y = 1;
    while y == 1 || weather.Timestamp(1)<=Date(end)
        index = Date<=weather.Timestamp(end) & Date>=weather.Timestamp(1);
        for j = 1:1:length(S)
            if isnumeric(weather.(S{j}))
                TestWeather.(S{j})(index,1) = interp1(weather.Timestamp,weather.(S{j}),Date(index)); 
            end
        end
        D3 = datevec(weather.Timestamp(1));
        weather.Timestamp = weather.Timestamp + datenum([D3(1)+1,1,1]) - datenum([D3(1),1,1]);%shift weather data 1 year
        y = y+1;
    end
    for j = 1:1:length(S)
        if isnumeric(weather.(S{j}))
            index = isnan(TestWeather.(S{j}));
            TestWeather.(S{j})(index,1) = weather.(S{j})(1);
        end
    end
end