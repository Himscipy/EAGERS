function [errors,Rsquare] = CompareForecastSimulate(building,weather,Date)
% Compare a forecasted building to it's simulated response
% building is the EAGERS building structure, weather is a structure of dry
% bulb, wetbuld and relative humididty, and solar irradiation,
% Date is the timsteps to be compared (must be sequential).

%% Handle input DATE
[s1, ~] = size(Date);
if s1 == 1
    Date = Date';
end
Location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
%% Run EAGERS method
Ti = 17.4;
iWeather = interpolateWeather(weather,Date);
Tdb = iWeather.Tdb;
RH = iWeather.RH;
[InternalGains,ExternalGains] = BuildingLoads(building,iWeather.irradDireNorm,iWeather.irradDiffHorz,Location,Date);

[Cooling, Heating, Fan_Power,Tzone,Twall,Damper] = BuildingProfile(building,Date,InternalGains,ExternalGains,Tdb,RH,Ti,Ti);


%% Simulate response to this Heating/Cooling Profile
AirFlow = Fan_Power/Build.VariableStruct.FanPower;
dt = Date(2:end) - Date(1:end-1);
dt = [dt(1);dt(1);dt];
[Tzone_S,Twall_S] = BuildingSimulate(building,Tdb,RH,dt,InternalGains,ExternalGains,Cooling, Heating,AirFlow,Damper,Tzone(1),Twall(1));

%% Plotting
errors = zeros(2, 1);
Rsquare = zeros(2, 1);
fprintf('Total %% errors:\n')



%compare air temp
figTitle = 'Zone Temperature';
fig = figure(1);
set(fig,'name',figTitle)
hold off
plotStep(Date,Tzone(2:end),'b')
hold on
plotStep(Date,Tzone_S(2:end),'r')
xlabel('Day of Year')
ylabel('Temperature (C)')
legend({'Forecast','Simulate'})
errors(1) = 100*(sum(Tzone_S.*dt)-sum(Tzone.*dt))/sum(Tzone.*dt);
fprintf('%s percent error:\t%f\n', figTitle, errors(1))
SSE = zeros(length(dt),1);
for i = 1:1:length(dt)
    SSE(i) = (sum(Tzone_S((i-1)+1:i).*dt((i-1)+1:i)) - Tzone(i).*dt(i))^2;
end
Rsquare(1) = 1 - sum(SSE)/sum((Tzone.*dt-mean(Tzone.*dt)).^2);
fprintf('%s R^2 error:\t%f\n', figTitle, Rsquare(1))

%compare interior lighting
figTitle = 'Wall Temperature';
fig = figure(2);
set(fig,'name',figTitle)
hold off
plotStep(Date,Twall(2:end),'b')
hold on
plotStep(Date,Twall_S(2:end),'r')
xlabel('Day of Year')
ylabel('Temperature (C)')
legend({'Forecast','Simulate'})
errors(2) = 100*(sum(Twall_S.*dt)-sum(Twall.*dt))/sum(Twall.*dt);
fprintf('%s percent error:\t%f\n', figTitle, errors(1))
SSE = zeros(length(dt),1);
for i = 1:1:length(dt)
    SSE(i) = (sum(Twall_S((i-1)+1:i).*dt((i-1)+1:i)) - Twall(i).*dt(i))^2;
end
Rsquare(2) = 1 - sum(SSE)/sum((Twall.*dt-mean(Twall.*dt)).^2);
fprintf('%s R^2 error:\t%f\n', figTitle, Rsquare(2))
end%Ends function CompareForecastSimulation

function plotStep(varargin)
Xe = varargin{1};
Ye = varargin{2};
if length(varargin) > 2
    lineSpec = {varargin{3:end}};
else
    lineSpec = {};
end
dt = Xe(2:end) - Xe(1:end-1);
[X,I] = sort([Xe(1)-dt(1);Xe;Xe(2:end)-dt+1e-9]);
Y = [Ye(1);Ye;Ye(2:end)];
Y = Y(I);
plot(X, Y, lineSpec{:})
end % Ends function plotStep
