function plotDemand(ForecastTime,Forecast,ActualTime,Actual)
%% plot demands, forecasted and actual into forecasting page of GUI
global Plant
handles = guihandles;
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
networkNames = fieldnames(Plant.subNet);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    networkNames(end+1) = {'BuildingTemp'};
end
nPlot = length(networkNames);
if isempty(ActualTime)
    Time = ForecastTime;
    Data = Forecast;
else
    Time = ActualTime;
    Data = Actual;
end
%% Make text strings for bottom axis
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Aug','Nov','Dec'};
D1 = datevec(Time(1));
D2 = datevec(Time(end));
nS = length(Time);
if floor(Time(1)) == floor(Time(round(2/3*nS)))
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))});
elseif floor(Time(1)) == floor(Time(end))-1% two days
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))},{'                       '},months(D2(2)),{' '},{num2str(D2(3))},{'  '},{num2str(D2(1))});
else %many days
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))},{'   through  '},months(D2(2)),{' '},{num2str(D2(3))},{'  '},{num2str(D2(1))});
end
D = datevec(Time(1));
hours = D(4) + D(5)/60 + D(6)/3600 + 24*(Time - Time(1));
if ~isempty(ActualTime) && ~isempty(ForecastTime)
    hours2 = D(4) + D(5)/60 + D(6)/3600 + 24*(ForecastTime - Time(1));
end
%% Do actual Plotting
for q = 1:1:nPlot
    h = handles.(strcat('ForecastPlot',num2str(q)));
    cla(h);
    if q==1
        tSize = 12;
    else
        tSize = 9;
    end
    S = get(handles.(strcat('ForecastName',num2str(q))),'String');
    if strcmp(S,'Electrical')
        S = 'E';
        ylab = 'Electrical Demand (kW)';
    elseif strcmp(S,'DistrictHeat')
        S = 'H';
        ylab = 'Heating Demand (kW)';
    elseif strcmp(S,'DistrictCool')
        S = 'C';
        ylab = 'Cooling Deamnd (kW)';
    elseif strcmp(S,'Hydro')
        S = 'W';
        ylab = 'Withdrawls (1000cfs)';
    elseif strcmp(S,'Steam')
        S = 'S';
        ylab = 'Steam Demand (1000 lbs/hr)';
    elseif strcmp(S,'BuildingTemp')
        S = 'T';
        ylab = 'Building Temperature Setpoint (C)';
    end
    if tSize==12 || isempty(ForecastTime)
        s0 = 1;
    else
        s0 = nnz(Time<=ForecastTime(1));
    end
    axTick = (ceil(hours(1)):round((hours(end)-hours(1))/12):hours(end));
    axIndex = mod(axTick,24);
    axIndex([false,axIndex(2:end)==0]) = 24;

    OoM = max(1,log10(max(sum(Data.(S)(s0:end,:),2))));
    if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
        Yspace = 10^(OoM-1);
        Ymax = 10^OoM;
    elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
        Yspace = 10^floor(OoM);
        Ymax = 10^ceil(OoM);
    elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
        Yspace = .5*10^floor(OoM);
        Ymax = .5*10^ceil(OoM);
    else  %count in increments of 2, 20, 200 or 2000 etc
        Yspace = .2*10^floor(OoM);
        Ymax = .2*10^ceil(OoM);
    end
    Ymin = 0;
    if strcmp(S,'T')
        Ymin = 55;
        Ymax = 85;
        Yspace = 5;
        Data.T = Data.T*9/5+32;
    end
    plot(h,hours(s0:end),Data.(S)(s0:end,:),'k')
    if ~isempty(ActualTime) && ~isempty(ForecastTime)
        plot(h,hours2,Forecast.(S),'b')
    end

    xlabel(h,dateText,'Color','k','FontSize',tSize)
    ylabel(h,ylab,'Color','k','FontSize',tSize)
    set(h,'XTick',axTick,'XTickLabel', {axIndex})
    if tSize==12 && ~isempty(ActualTime) && ~isempty(ForecastTime)
        a = hours(nnz(ActualTime<=ForecastTime(1)));
        if a>1
            plot(h,[a,a],[Ymin,Ymax],'c--')
        end
    end
    xlim(h,[hours(s0), hours(end)])
    ylim(h,[Ymin,Ymax])
    set(h,'YTick',Ymin:Yspace:Ymax,'FontSize',tSize-2)
    h2 = handles.(strcat('ForecastPlot',num2str(q),'b'));
    set(h2,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel', [])
    ylabel(h2,[]);
end