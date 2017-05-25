function plotDispatch(handles,ForecastTime,Forecast,HistoryTime,History)
%% plot dispatch into GUI, with historical operation
global Plant
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
nG = length(Plant.Generator);
stor = [];
Names = cell(nG,1);
for i = 1:1:nG
    Names(i) = {Plant.Generator(i).Name};
    if ~isempty(strfind(Plant.Generator(i).Type,'Storage')) && Plant.Generator(i).Enabled
        stor(end+1) = i;
    end
end
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
nPlot = length(networkNames);

%% Collate history and future data  
if isempty(History)
    backSteps = 0;
    Data = Forecast;
    Time = ForecastTime*24;
else
    backSteps = length(History(:,1));
    Data = [History;Forecast];
    Time = [HistoryTime;ForecastTime]*24;
end
D = datevec(Time(1)/24);
hours = D(4)+D(5)/60+D(6)/3600+(Time - Time(1));
dt = Time(2:end) - Time(1:end-1);

%% Make text strings to scroll across bottom axis
Dprev = datevec(ForecastTime(1)-1);
Dnext = datevec(ForecastTime(1)+1);
Dtimestamp = ForecastTime(1)+hours/24;
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Aug','Nov','Dec'};
dateText = ' ';
if backSteps<24/Plant.optimoptions.Resolution
    b = ceil(30*(24/Plant.optimoptions.Resolution-backSteps)/(24/Plant.optimoptions.Resolution));
    for i = 1:1:b
        dateText = strcat(dateText,{' '});
    end
end
if nnz(Dtimestamp<datenum([D(1),D(2),D(3)]))>0.3*length(hours) %more than 30% of window is previous day
    dateText = strcat(dateText,months(Dprev(2)),{' '},{num2str(Dprev(3))},{'  '},{num2str(Dprev(1))},{'                                '});
else  %append spaces
    c = ceil(30*nnz(Dtimestamp<datenum([D(1),D(2),D(3)]))/length(hours));
    for i = 1:1:c
        dateText = strcat(dateText,{' '});
    end
end
dateText = strcat(dateText,months(D(2)),{' '},{num2str(D(3))},{'  '},{num2str(D(1))});
if nnz(Dtimestamp>datenum([Dnext(1),Dnext(2),Dnext(3)]))>0.3*length(hours) %more than 30% of window is next day
    dateText = strcat(dateText,{'                         '},months(Dnext(2)),{' '},{num2str(Dnext(3))},{'  '},{num2str(Dnext(1))});
else  %append spaces
        dateText = strcat(dateText,{'                                      '});
end
dateText2 = strcat('  ',months(D(2)),{' '},{num2str(D(3))},{'  '},{num2str(D(1))});
if nnz(Dtimestamp>datenum([Dnext(1),Dnext(2),Dnext(3)]))>0.3*length(hours) %more than 30% of window is next day
    dateText2 = strcat(dateText2,{'                      '},months(Dnext(2)),{' '},{num2str(Dnext(3))},{'  '},{num2str(Dnext(1))});
else
    dateText2 = strcat(dateText2,{'                                '});
end

%% convert the saved SOC to power
StoragePower = 0*Data;
StorageState = 0*Data;
for i = 1:1:length(stor)
    StorageState(:,stor(i)) = Data(:,stor(i))+ ones(length(hours),1)*(Plant.Generator(stor(i)).OpMatA.Stor.Size - Plant.Generator(stor(i)).OpMatA.Stor.UsableSize); %add the unusable charge
    if strcmp(Plant.Generator(stor(i)).Type,'Hydro Storage')%% Need the river segment and spill flow to calculate power
        StoragePower(2:end,stor(i)) = (Data(2:end,Plant.Generator(stor(i)).OpMatA.DownRiverSegment) - Data(2:end,Plant.Generator(stor(i)).OpMatA.SpillFlow))*Plant.Generator(stor(i)).OpMatA.output.E;
    else
        StoragePower(2:end,stor(i)) = (StorageState(1:end-1,stor(i)) - StorageState(2:end,stor(i)))./dt;  
    end
end
Data(:,stor) = StoragePower(:,stor);

%% Do actual Plotting
colorVec = get(handles.figure1,'UserData');
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,length(Names)));
if get(handles.StackedGraph,'Value')==0 || strcmp(get(handles.uipanelMain3,'visible'),'on')
    LINE = true;
else
    LINE = false;
end
for q = 1:1:nPlot
    if strcmp(get(handles.(strcat('ResultPlot',num2str(q))),'Visible'),'on')
        h = handles.(strcat('ResultPlot',num2str(q)));
        cla(h);
        h2 = handles.(strcat('ResultPlot',num2str(q),'b'));
        cla(h2);
        s1 = max(1,length(HistoryTime));
        if q==1
            tSize = 12;
        else
            tSize = 9;
        end
        S = get(handles.(strcat('ResultName',num2str(q))),'String');
        if strcmp(S,'Electrical')
            S = 'E';
        elseif strcmp(S,'DistrictHeat')
            S = 'H';
        elseif strcmp(S,'DistrictCool')
            S = 'C';
        elseif strcmp(S,'Hydro')
            S = 'W';
            for i = 1:1:length(Plant.Generator)
                if strcmp(Plant.Generator(i).Type,'Hydro Storage')%% plot the downriver flow
                    Data(:,i) = Data(:,Plant.Generator(i).OpMatA.DownRiverSegment);
                end
            end
        elseif strcmp(S,'Steam')
            S = 'S';
        end

        [name,stor,Data2,storNames] = sortForPlot(S,Data);
        %% Storage
        if length(stor)>2
            if strcmp(get(handles.(strcat('ResultStorageSelect',num2str(q))),'Visible'),'off')
                set(handles.(strcat('ResultStorageSelect',num2str(q))),'Visible','on','value',1,'string',storNames)
            end
            storNum = get(handles.(strcat('ResultStorageSelect',num2str(q))),'value');%select only 1 storage at a time
            storRmv = storNames(storNum~=stor);
            stor = stor(storNum); 
            storNames = storNames(storNum);
            I = 1:length(name);
            for r = 1:1:length(name)
                if any(strcmp(name(r),storRmv))
                    I(r) = 0;
                end
            end
            I = nonzeros(I);
            name = name(I);
            Data2 = Data2(:,I);
        else
            if strcmp(get(handles.(strcat('ResultStorageSelect',num2str(q))),'Visible'),'on')
                set(handles.(strcat('ResultStorageSelect',num2str(q))),'Visible','off')
            end
        end
        %% Plot    
        if ~isempty(Data2)
            dataPlot(h,LINE,Data2,hours,dt,s1,name,Names,colorsPlot,tSize)
            
            if ~isempty(stor) && ~LINE
                addStorage2Plot(h,h2,StorageState(:,stor),hours,s1,colorsPlot,storNames,Names,tSize,handles.(strcat('ResultName',num2str(q))))
            else
                set(h2,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[])
                ylabel(h2,[]);
            end
            if q==1
                xlabel(h,dateText,'Color','k','FontSize',tSize)
            else xlabel(h,dateText2,'Color','k','FontSize',tSize)
            end
            if strcmp(get(handles.(strcat('ResultName',num2str(q))),'String'),'Hydro')
                ylabel(h,'Flow (1000cfs)','Color','k','FontSize',tSize)
            else
                ylabel(h,'Generation (kW)','Color','k','FontSize',tSize)
            end
        end
    end
end
A = {ForecastTime,Forecast,HistoryTime,History};
set(handles.ResultPlot1,'UserData',A);

function [name,stor,Data2,storNames] = sortForPlot(S,Data)
global Plant
name = {};
storNames ={};
stor =[];
Data2 =[];
for i = 1:1:length(Plant.Generator)
    include = false;
    dataI = Data(:,i);
    if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S) && Plant.Generator(i).Enabled % energy storage
        if strcmp(S,'E') && strcmp(Plant.Generator(i).Type,'Hydro Storage')
            %don't plot state-of charge on electric plot
        else
            stor(end+1) = i;
        end
        include = true;
    elseif (strcmp(Plant.Generator(i).Type,'Utility')||~isempty(strfind(Plant.Generator(i).Type,'District'))) && isfield(Plant.Generator(i).OpMatA.output,S) %utilities
        include = true;
    elseif isfield(Plant.Generator(i).OpMatA.output,S) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S,'E')) %generators
        include = true;
        if strcmp(S,'H') && isfield(Plant.Generator(i).OpMatA.output,'E')
            dataI = dataI*Plant.Generator(i).OpMatA.output.(S)(1); %Hratio for CHP generators
        end
    elseif strcmp(S,'E') && strcmp(Plant.Generator(i).Source,'Renewable')
        include = true;
    end
    if include
        name(end+1) = {Plant.Generator(i).Name};
        Data2(:,end+1) = dataI;
        if any(stor==i)% energy storage
            storNames(end+1) ={Plant.Generator(i).Name};
        end
    end
end

function dataPlot(h,LINE,Data2,hours,dt,s1,name,Names,colorsPlot,tSize)
if tSize==12
    s0 = 1;
else s0 = s1;
end
axTick = (ceil(hours(1)):round((hours(end)-hours(1))/12):hours(end));
axIndex = mod(axTick,24);
axIndex([false,axIndex(2:end)==0]) = 24;

plotTime = zeros(2*length(hours(s0:end))-2,1);
plotTime(1:2:2*length(hours(s0:end))-2) = [hours(s0);hours(s0+2:end)-.9999*dt(s0+1:end)];
plotTime(2:2:2*length(hours(s0:end))-2) = hours(s0+1:end);
Data3 = zeros(length(Data2(s0:end,1))*2-2,length(Data2(1,:)));
Data3(1:2:end,:) = Data2(s0+1:end,:);
Data3(2:2:end,:) = Data2(s0+1:end,:);
if LINE
    str2 = 'Color';
    h1 = plot(h,plotTime,Data3,'LineWidth',3);
    negBars = [];
else
    str2 = 'FaceColor';
    posBars = max(0,Data3);
    negBars = min(0,Data3);
    Ikeep = any(negBars<0);
    negBars = negBars(:,Ikeep);
    negBarNames = name(Ikeep);
    h1 = area(h,plotTime,posBars,'Linestyle','none');
    if ~isempty(negBars)
        L = area(h,plotTime,negBars,'Linestyle','none');
        for c = 1:1:length(L)
            set(L(c),str2,colorsPlot(strcmp(negBarNames(c),Names),:));
        end
    end
end
for c = 1:1:length(h1)
    set(h1(c),str2,colorsPlot(strcmp(name(c),Names),:));
end

if ~isempty(negBars)
    OoM = log10(max(sum(posBars,2)-sum(negBars,2)));
else OoM = log10(max(sum(Data3,2)));
end
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

negTicks = floor(min(min(min(0,Data3)))/Yspace);
if abs(negTicks)>3
    Yspace = 2*Yspace;
    negTicks = floor(min(min(min(0,Data3)))/Yspace);
end
Ymin = Yspace*negTicks;
if isempty(Ymin)
    Ymin = 0;
end
if tSize==12 && s1>1
    plot(h,[hours(s1),hours(s1)],[Ymin,Ymax],'c--')   
end
xlim(h,[hours(s0), hours(end)])
ylim(h,[Ymin,Ymax])
axTickY = Ymin:Yspace:Ymax;
set(h,'YTick',axTickY,'YTickLabel', {axTickY},'FontSize',tSize-2)
set(h,'XTick',axTick,'XTickLabel', {axIndex})

function addStorage2Plot(h,h2,StorageState,hours,s1,colorsPlot,storNames,Names,tSize,hName)
if tSize==12
    s0 = 1;
else s0 = s1;
end

L = plot(h2,hours(s0:end),StorageState(s0:end,:));
for c = 1:1:length(L)
    set(L(c),'Color',colorsPlot(strcmp(storNames(c),Names),:),'LineStyle','-','LineWidth',2,'Marker','x','MarkerEdgeColor','k','MarkerSize',5)
end
xlim(h2,get(h,'Xlim'));
set(h2,'xtick',[],'xticklabel',[])
ticks = get(h,'YTick');

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
negTicks = min(ticks)/(ticks(end)-ticks(end-1));
Ymin = Yspace*negTicks;  
ylim(h2,[Ymin,Ymax])
axTickY = Ymin:Yspace:Ymax;
set(h2,'YTick',axTickY,'YTickLabel', {axTickY},'FontSize',tSize-2)
if strcmp(get(hName,'String'),'Hydro')
    ylabel(h2,'State of Charge (1000 acre-ft)','Color','k','FontSize',tSize)
else
    ylabel(h2,'State of Charge (kWh)','Color','k','FontSize',tSize)
end
if ~isempty(storNames)
    legend(h2,storNames,'Fontsize',tSize-1,'Orientation','Horizontal','Location','North','Box','off','color','none');%,'Xcolor',[1 1 1],'Ycolor',[1 1 1])
end