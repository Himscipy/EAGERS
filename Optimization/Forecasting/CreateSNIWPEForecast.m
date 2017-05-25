function Forecast = CreateSNIWPEForecast(Date,RealData)
global Last24hour

%% initialize
beta = 0.9; % Take this out after training is added!
S = fieldnames(Last24hour.Demand);
Forecast.Timestamp = Date(2:end); % first value of Date is current time

%% forecast horizons and historical data times
hrstep = 1/24; % equivalent to 1 hour in datenum format
nmlzDate = Date - Date(1);
hrDate = nmlzDate ./ hrstep;
hzn = hrDate(2:end); % first value of Date is current time

histNum = [Last24hour.Timestamp;Date(1)]; % include current datum
nmlzhistNum = histNum - histNum(1);
histTime = nmlzhistNum ./ hrstep;

%% check for existence of prevErr and prevFcast in Last24hour
if ~isfield(Last24hour,'prevErr')
    Last24hour.prevErr.Temp = 0;
    for i = 1:1:length(S)
        Last24hour.prevErr.(S{i}) = 0;
    end
end
if ~isfield(Last24hour,'prevFcast')
    Last24hour.prevFcast.Temp = RealData.Temperature;
    for i = 1:1:length(S)
        Last24hour.prevFcast.(S{i}) = RealData.Demand.(S{i});
    end
end

%% create forward forecast of temperature
histData = [Last24hour.Temperature;RealData.Temperature];
pf = Last24hour.prevFcast.Temp;
pe = Last24hour.prevErr.Temp;
[Forecast.Temperature, Last24hour.prevErr.Temp] = ...
    SNIWPEForecast(histData, histTime, hzn, pf, pe, beta);
Last24hour.prevFcast.Temp = Forecast.Temperature(1);

%% repeat for electric, cooling, heating, and steam demands as necessary
for i = 1:1:length(S)
    histData = [Last24hour.Demand.(S{i});RealData.Demand.(S{i})];
    pf = Last24hour.prevFcast.(S{i});
    pe = Last24hour.prevErr.(S{i});
    [Forecast.Demand.(S{i}), Last24hour.prevErr.(S{i})] = ...
        SNIWPEForecast(histData, histTime, hzn, pf, pe, beta);
    Last24hour.prevFcast.(S{i}) = Forecast.Demand.(S{i})(1);
end
