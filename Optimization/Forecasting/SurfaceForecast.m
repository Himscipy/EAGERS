function Forecast = SurfaceForecast(Date,Tforecast)
global Plant Last24hour
dSim = Last24hour.Timestamp(end);
Days = ceil(Date(end)-floor(dSim)); %this says that you may need to use surface fits for multiple days.
s12hour = nnz(Date<(dSim+.5)); %steps in next 12 hours
A = datevec(dSim);
dateDays = floor(dSim)-1:ceil(dSim+Days);
Weekday = zeros(length(dateDays),1);
for i = 1:1:length(dateDays)
    % Is dateDays(i) a weekday and not a holiday?
    if weekday(dateDays(i))<=6 && weekday(dateDays(i))>=2 && (~isfield(Plant.Data,'Holidays') || ~any(Plant.Data.Holidays==dateDays(i)))
        Weekday(i) = 1;
    else
        Weekday(i) = 0;
    end
    % Is the historical data from a year in the past? (DOES NOT HANDLE
    % HISTORICAL DATA FROM MORE THAN 1 YEAR IN THE PAST.)
    if Weekday(i)==1 && (~isfield(Plant.Data,'Holidays') || (any(Plant.Data.Holidays == dateDays(i)+365) && dateDays(i)<min(Plant.Data.Holidays)))
        Weekday(i) = 0;
    end
end
monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};
hour = (Date-floor(dSim))*24;
S = fieldnames(Last24hour.Demand);
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [~,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.Demand.(S{s}) = zeros(length(Date),nD);
    list = fieldnames(Plant.Data.HistProf.(S{s})(1));
    for k = 1:1:nD
        month = A(2);
        HistFitDem = nan(length(Date),1);
        for j = 0:1:ceil(Days)
            if A(3)+ (j-1) > (datenum(A(1), A(2)+1, 1)-datenum(A(1), A(2), 1))
                month = A(2)+1;
            end
            %load suface
            if length(list)>1
                if nnz(strcmp(list,strcat(monthName(month),'WeekEnd')))>0 %historical profile is split to weekday/weekend
                    if Weekday(j+1) == 0
                        Surface = Plant.Data.HistProf.(S{s})(k).(strcat(monthName{month},'WeekEnd'));
                    else
                        Surface = Plant.Data.HistProf.(S{s})(k).(strcat(monthName{month},'WeekDay'));
                    end
                elseif nnz(strcmp(list,char(monthName(month))))>0 %historical profile is broken by month, but not weekend/weekday
                    Surface = Plant.Data.HistProf.(S{s})(k).(monthName{month});
                else 
                    Surface = Plant.Data.HistProf.(S{s})(k).(list{1});
                end
            else
                Surface = Plant.Data.HistProf.(S{s})(k).(list{1});% only 1 historical surface fit
            end
            if j ==0 %surface for yesterday
                YestPred = Surface(mod((Last24hour.Timestamp-floor(Last24hour.Timestamp(1)))*24,24),Last24hour.Weather.Tdb);
                errorPred = min(1,max(-1,(Last24hour.Demand.(S{s})(:,k)- YestPred)./YestPred));
                BiasPerc = mean(errorPred)*100; %average percent error between hisorical expectations for yesterday, and actual yesterday
                HourlyError = Last24hour.Demand.(S{s})(:,k) - YestPred*(1+0.7*BiasPerc/100);%hourly variation with equal mean (historical and last 24 hours)
            else %surface for current day or subsequent days
                index = (hour<=(24*j)) & (hour>=(24*(j-1))); %index of Time corresponding to this day
                HistFitDem(index) = Surface(mod(hour(index),24),Tforecast(index));
            end
        end
        Forecast.Demand.(S{s})(1:s12hour,k) = max(0,HistFitDem(1:s12hour)*(1+0.7*BiasPerc/100) + linspace(1,0,s12hour)'.*interp1(Last24hour.Timestamp,HourlyError,Date(1:s12hour)-1));
        Forecast.Demand.(S{s})(s12hour+1:end,k) = max(0,HistFitDem(s12hour+1:end)*(1+0.7*BiasPerc/100));
        Forecast.Demand.(S{s})(abs(Forecast.Demand.(S{s}))<1e-2) = 0; 
%         if strcmp(Plant.Name,'PNNL_C') && strcmp(S{s},'C') && any(Forecast.Demand.C<0)
% %             disp('Forecast is negative')
%             Forecast.Demand.C = max(0,Forecast.Demand.C);
% %             Forecast.Demand.(S{s})(:,k) = interp1(Last24hour.Timestamp,Last24hour.Demand.(S{s})(:,k),Date-1);%for chilling overwrite with yesterdays cooling demand
%         end
    end
end