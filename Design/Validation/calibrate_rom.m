function [errors,Rsquare] = calibrate_rom(building,weather,Date)
% Compare a forecasted building to it's energyPlus counterpart
% building is the EAGERS building structure, weather is a structure of dry
% bulb, wetbuld and relative humididty, and solar irradiation,
% Date is the timsteps to be compared (must be sequential).

%% Handle input DATE
[s1, ~] = size(Date);
if s1 == 1
    Date = Date';
end
Location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
%% Load EPlus meter file
global Model_dir
load(fullfile(Model_dir,'Design','Validation',strcat(building.Name,'.mat')))% building
load('Tbuild.mat');%temperature from EnergyPlus
iWeather = interpolate_weather(weather,Date);
Tdb = iWeather.Tdb;
RH = iWeather.RH;
D = datevec(Date(1));
DofY = Date - datenum([D(1), 1,1]);
dt = Date(2:end) - Date(1:end-1);
dt = [dt(1);dt];
if isempty(mtr.Timestamp)
    EplusTime = linspace(datenum([2017,1,1,1,0,0]),datenum([2018,1,1]),8760)';
else
    EplusTime = datenum(datevec(mtr.Timestamp));
end
D = datevec(EplusTime(1));
EplusTime = EplusTime - datenum([D(1), 1,1]);
xi = max(1,nnz((EplusTime)<=DofY(1)));
xf = nnz((EplusTime)<=DofY(end));

SG = solar_gain(building,Date,Location,iWeather);
ExternalGains = SG.Walls + SG.Roof;
B_loads = building_loads(building,Date,SG);
Cooling = zeros(length(Tdb),1);
Heating = zeros(length(Tdb),1);
AirFlow = zeros(length(Tdb),1);
Tzone = zeros(length(Tdb)+1,1);
Twall = zeros(length(Tdb)+1,1);
for i = 1:1:(xf-xi)+1
    Cooling(1+10*(i-1):10*i) = mtr.CoolingElectricity(i+xi-1)*building.VariableStruct.COP_C;
    Heating(1+10*(i-1):10*i) = mtr.HeatingElectricity(i+xi-1)*building.VariableStruct.COP_H + mtr.HeatingGas(i+xi-1);
    AirFlow(1+10*(i-1):10*i) = mtr.FansElectricity(i+xi-1)/building.VariableStruct.FanPower;
end
Damper = 0.1*ones(length(AirFlow),1);
Tzone(1) = 17;
Twall(1) = 17;
for t = 1:1:length(Tdb)
    [zone,wall] = building_simulate(building,Tdb(t),RH(t),dt(t)*24*3600,B_loads.InternalGains(t),ExternalGains(t),Cooling(t),Heating(t),AirFlow(t),Damper(t),Tzone(t),Twall(t));
    if rem(t,240) == 0
        Tzone(t+1) = Tbuild(t);
        Twall(t+1) = Tbuild(t);
    else
        Tzone(t+1) = zone(end);
        Twall(t+1) = wall(end);
    end
end

%% Plotting
errors = zeros(1, 1);
Rsquare = zeros(1, 1);
fprintf('Total %% errors:\n')

%compare zone temp

figTitle = 'Zone Temperature';
ylab = 'Temperature (C)';
figNum = 1;
ePlus = Tbuild;
eagers = Tzone(2:end);
[errors(figNum),Rsquare(figNum)] = elec_load_plot_setup(figTitle, figNum, DofY, ePlus, DofY, eagers, dt, dt,ylab);

end % Ends function calibrate_rom


function [PercError,Rsquare] = elec_load_plot_setup(figTitle, figNum, Xe, E_Plus,DofY, Eagers, dt, dtE,ylab)
fig = figure(figNum);
set(fig,'name',figTitle)
hold off
plotStep(Xe,E_Plus,'b')
hold on
plotStep(DofY,Eagers,'r')
xlabel('Day of Year')
ylabel(ylab)
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Eagers.*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
fprintf('%s percent error:\t%f\n', figTitle, PercError)
n = round(dtE(1)/dt(1));
SSE = zeros(length(dtE),1);
for i = 1:1:length(dtE)
    SSE(i) = (sum(Eagers((i-1)*n+1:i*n).*dt((i-1)*n+1:i*n)) - E_Plus(i).*dtE(i))^2;
end
Rsquare = 1 - sum(SSE)/sum((E_Plus.*dtE-mean(E_Plus.*dtE)).^2);
fprintf('%s R^2 error:\t%f\n', figTitle, Rsquare)
end % Ends function elec_load_plot_setup


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
