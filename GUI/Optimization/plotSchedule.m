function plotSchedule(ForecastTime,Forecast,HistoryTime,History,ActualTime,Actual,gen)
%% plot dispatch into GUI, with historical operation
global Plant
handles = guihandles;
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
stor = false;

Name = Plant.Generator(gen).Name;
if ~isempty(strfind(Plant.Generator(gen).Type,'Storage'))
    stor = true;
end
Outs = fieldnames(Plant.Generator(gen).OpMatA.output);

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));
networkNames = networkNames(~strcmp('Buildings',networkNames));
nPlot = length(networkNames);

%% Collate history and future data
if isempty(History)
    Data = Forecast;
    Time = ForecastTime*24;
else
    Data = [History;Forecast];
    Time = [HistoryTime;ForecastTime]*24;
end
dt = Time(2:end) - Time(1:end-1);

%% Make text strings for bottom axis
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Aug','Nov','Dec'};
if ~isempty(ActualTime)
    Date = ActualTime;
else Date = [HistoryTime;ForecastTime;];
end
D1 = datevec(Date(1));
if floor(Date(1)) == floor(Date(end))
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))});
elseif floor(Date(1)) == floor(Date(end))-1% two days
    D2 = datevec(Date(end));
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))},{'                       '},months(D2(2)),{' '},{num2str(D2(3))},{'  '},{num2str(D2(1))});
else %many days
    D2 = datevec(Date(floor(length(Date)/2)));
    D3 = datevec(Date(end));
    dateText = strcat(months(D1(2)),{' '},{num2str(D1(3))},{'  '},{num2str(D1(1))},{'         '},months(D2(2)),{' '},{num2str(D2(3))},{'  '},{num2str(D2(1))},{'            '},months(D3(2)),{' '},{num2str(D3(3))},{'  '},{num2str(D3(1))});
end

D = datevec(Date(1));
hours = D(4) + D(5)/60 + D(6)/3600 + 24*(Date - Date(1));
hours2 = D(4) + D(5)/60 + D(6)/3600 + Time - 24*Date(1);
%% convert the saved SOC to power
StoragePower = 0*Data;
if stor
    StorageState = Data + ones(length(hours),1)*(Plant.Generator(gen).OpMatA.Stor.Size - Plant.Generator(gen).OpMatA.Stor.UsableSize); %add the unusable charge
    StoragePower(2:end) = (StorageState(1:end-1) - StorageState(2:end))./dt;
    Data = StoragePower;
end


%% Do actual Plotting
for q = 1:1:nPlot
    h = handles.(strcat('ForecastPlot',num2str(q)));
    cla(h);
    h2 = handles.(strcat('ForecastPlot',num2str(q),'b'));
    cla(h2);
    if q==1
        tSize = 12;
    else
        tSize = 9;
    end
    S = get(handles.(strcat('ForecastName',num2str(q))),'String');
    if strcmp(S,'Electrical')
        S = 'E';
    elseif strcmp(S,'DistrictHeat')
        S = 'H';
    elseif strcmp(S,'DistrictCool')
        S = 'C';
    elseif strcmp(S,'Hydro')
        S = 'W';
    elseif strcmp(S,'Steam')
        S = 'S';
    end
    if strcmp(Outs(1),S)
        set(h,'Visible','on');
        set(h2,'Visible','on');
        set(handles.(strcat('ForecastName',num2str(q))),'Visible','on');
        if tSize==12
            s0 = 1;
        else s0 = max(1,nnz(ActualTime<ForecastTime(1)));
        end
        axTick = (ceil(hours(1)):round((hours(end)-hours(1))/12):hours(end));
        axIndex = mod(axTick,24);
        axIndex([false,axIndex(2:end)==0]) = 24;

        OoM = log10(max(Forecast(s0:end)));
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
        negTicks = floor(min(min(0,Data))/Yspace);
        if abs(negTicks)>3
            Yspace = 2*Yspace;
            negTicks = floor(min(min(0,Data))/Yspace);
        end
        Ymin = Yspace*negTicks;
        if isempty(Ymin)
            Ymin = 0;
        end
        if ~isempty(Actual)
            plot(h,hours(s0:end),Actual(s0:end),'k')
        end
        plot(h,hours2,Data,'b')

        xlabel(h,dateText,'Color','k','FontSize',tSize)
        if strcmp(S,'W')
            ylabel(h,'Flow (1000cfs)','Color','k','FontSize',tSize)
        else
            ylabel(h,'Generation (kW)','Color','k','FontSize',tSize)
        end
        set(h,'XTick',axTick,'XTickLabel', {axIndex})
        xlim(h,[hours(s0), hours(end)])
        ylim(h,[Ymin,Ymax])
        set(h,'YTick',Ymin:Yspace:Ymax,'FontSize',tSize-2)
        if ~stor
            legend(h,{Name},'Fontsize',tSize-1,'Orientation','Horizontal','Location','Best','Box','off','color','none');%,'Xcolor',[1 1 1],'Ycolor',[1 1 1])
            set(h2,'YTick',[],'YTickLabel', [])
            ylabel(h2,[]);
        else
            legend(h,{[Name,'  Output']},'Fontsize',tSize-1,'Orientation','Horizontal','Location','Best','Box','off','color','none');%,'Xcolor',[1 1 1],'Ycolor',[1 1 1])
            L =plot(h2,hours2,StorageState);
            set(L,'Color','red','LineStyle','-','LineWidth',2,'Marker','x','MarkerEdgeColor','k','MarkerSize',5)
            xlim(h2,[hours(s0), hours(end)])
            OoM = log10(max(max(StorageState)));
            if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
                Ymax = 10^OoM;
            elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
                Ymax = 10^ceil(OoM);
            elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
                Ymax = .5*10^ceil(OoM);
            else  %count in increments of 2, 20, 200 or 2000 etc
                Ymax = .2*10^ceil(OoM);
            end
            pTicks = nnz(get(h,'YTick')>0); % # of positive tick marks
            Yspace = Ymax/pTicks;
            negTicks = Ymin/Yspace;
            Ymin = Yspace*negTicks;  
            ylim(h2,[Ymin,Ymax])
            axTickY = Ymin:Yspace:Ymax;
            set(h2,'xtick',[],'xticklabel',[],'YTick',axTickY ,'YTickLabel', {axTickY})
            ylabel(h2,'State of Charge (kWh)','Color','k','FontSize',tSize)
            legend(h2,{'State of Charge'},'Fontsize',tSize-1,'Orientation','Horizontal','Location','Best','Box','off','color','none');%,'Xcolor',[1 1 1],'Ycolor',[1 1 1])
        end
    else
        set(h,'Visible','off');
        set(h2,'Visible','off');
        set(handles.(strcat('ForecastName',num2str(q))),'Visible','off');
    end
end