function plotDispatch(handles,ForecastTime,Solution,HistoryTime,History)
%% plot dispatch into GUI, with historical operation
global Plant CurrentState
if isfield(handles,'LegendDeleteProxy')%2013 matlab
%     delete(handles.LegendColorbarLayout)
%     delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
stor = [];
Names = cell(nG+nB,1);
for i = 1:1:nG
    Names(i) = {Plant.Generator(i).Name};
    if isfield(Plant.Generator(i).QPform,'Stor') && Plant.Generator(i).Enabled
        stor(end+1) = i;
    end
end
for i = 1:1:nB
    Names(nG+i) = {Plant.Building(i).Name};
end

%% Collate history and future data
if isempty(History.Dispatch)
    backSteps = 0;
    Data = Solution.Dispatch;
    Time = ForecastTime*24;
else
    backSteps = length(History.Dispatch(:,1));
    if abs(ForecastTime(1) - HistoryTime(end))<1e-6
        Data = [History.Dispatch;Solution.Dispatch(2:end,:)];
        Time = [HistoryTime;ForecastTime(2:end)]*24;
    else
        Data = [History.Dispatch;Solution.Dispatch];
        Time = [HistoryTime;ForecastTime]*24;
    end
end
D = datevec(Time(1)/24);
hours = D(4)+D(5)/60+D(6)/3600+(Time - Time(1));
dt = Time(2:end) - Time(1:end-1);

%% Make text strings to scroll across bottom axis
Dprev = datevec(ForecastTime(1)-1);
Dnext = datevec(ForecastTime(1)+1);
Dtimestamp = ForecastTime(1)+hours/24;
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
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
dateText = {dateText;dateText2;};
%% convert the saved SOC to power
StoragePower = 0*Data;
StorageState = 0*Data;
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        StoragePower(:,i) =  Data(:,i);
        StorageState(1,i) = CurrentState.Hydro(1,i);
        if isempty(History.hydroSOC)
            StorageState(2:end,i) = Solution.hydroSOC(:,i);
        else
            StorageState(:,i) = [History.hydroSOC(:,i);Solution.hydroSOC(:,i);];
        end
    elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydrogen Storage'})
        StorageState(:,i) = Data(:,i);
        if isfield(Plant.Generator(i).VariableStruct,'MaxDOD')
            StorageState(:,i) = Data(:,i)+ Plant.Generator(i).Size*(1-Plant.Generator(i).VariableStruct.MaxDOD/100); %add the unusable charge
        end
        StoragePower(2:end,i) = (StorageState(1:end-1,i) - StorageState(2:end,i))./dt;  
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
nPlot = 0;
while isfield(handles,strcat('ResultPlot',num2str(nPlot+1)))
    nPlot = nPlot+1;
    if strcmp(get(handles.(strcat('ResultPlot',num2str(nPlot))),'Visible'),'on')
        h = handles.(strcat('ResultPlot',num2str(nPlot)));
        cla(h);
        h2 = handles.(strcat('ResultPlot',num2str(nPlot),'b'));
        cla(h2);
        s1 = max(1,length(HistoryTime));
        if nPlot==1
            tSize = 12;
        else
            tSize = 9;
        end
        S = get(handles.(strcat('ResultName',num2str(nPlot))),'String');
        if strcmp(S,'BuildingTemp')
            name = cell(nB,1);
            storNames ={};
            stor =[];
            for i = 1:1:nB
                name(i) = {Plant.Building(i).Name};
            end
            if isempty(History.Buildings)
                Data = [CurrentState.Buildings(1,:);Solution.Buildings.Temperature]*9/5+32;
            else
                Data = [History.Buildings;Solution.Buildings.Temperature]*9/5+32;
            end
            dataPlot(h,nPlot,true,Data,hours,dt,s1,name,Names,colorsPlot,tSize,S,dateText)
            ylim(h,[55,80])
            set(h2,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[])
            ylabel(h2,[]);
            xlabel(h,dateText2,'Color','k','FontSize',tSize)
            ylabel(h,'Building Temperature Setpoint (F)','Color','k','FontSize',tSize)
            Data2 = [];
        elseif strcmp(S,'Hydro')
            stor = [];
            storNames ={};
            for i = 1:1:length(Plant.Generator)
                if strcmp(Plant.Generator(i).Type,'Hydro Storage')%% plot the downriver flow
                    n = Plant.Generator(i).QPform.DownRiverSegment;
                    if isempty(History.LineFlows)
                        Data2(1,i) = 0;
                        Data2(2:end,i) = Solution.LineFlows(:,n);%Outflow is after Power Produced
                    else
                        Data2(:,i) = [History.LineFlows(:,n);Solution.LineFlows(:,n)];%Outflow is after Power Produced
                    end
                    stor(end+1) = i;
                    storNames(end+1) = {Plant.Generator(i).Name};
                end
            end
        else
            if ~isfield(Plant.Generator(i),'QPform')%NREL case
                name = {'Elec Utility1';'fuel cell';'Battery';};
                stor = 5;
                Data2 = [Data(:,1),Data(:,3),Data(:,5)];
                storNames = {'Battery'};
            else 
                [name,stor,Data2,storNames] = sortForPlot(S,Data);
            end
        end
        %% Plot    
        if ~isempty(Data2)
            if ~isempty(stor) && ~LINE
                PlotWithStorage(Data2,StorageState(:,stor),dt,hours(end-length(StorageState(:,1))+1:end,:),s1,colorsPlot,storNames,Names,tSize,stor,name,nPlot,handles,S,dateText)
            else
                dataPlot(h,nPlot,LINE,Data2,hours,dt,s1,name,Names,colorsPlot,tSize,S,dateText)
                set(h2,'xtick',[],'xticklabel',[],'YTick',[],'YTickLabel',[])
                ylabel(h2,[]);
            end
        end
    end
end
Solution.ForecastTime = ForecastTime;
Solution.HistoryTime = HistoryTime;
Solution.History = History;
set(handles.ResultPlot1,'UserData',Solution);
end%Ends function plotDispatch

function [name,stor,Data2,storNames] = sortForPlot(S,Data)
global Plant
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';};
nnAbrev = {'E';'H';'C';'W';'DC';'CW';'E1';'Hy';'LH2';'H2';};
index = nonzeros((1:10)'.*strcmp(S,nnList));
S = nnAbrev{index};
name = {};
storNames ={};
stor =[];
Data2 =[];
for i = 1:1:length(Plant.Generator)
    include = false;
    if isfield(Plant.Generator(i).QPform,'Stor') && isfield(Plant.Generator(i).QPform.output,S) && Plant.Generator(i).Enabled % energy storage
        if strcmp(S,'E') && strcmp(Plant.Generator(i).Type,'Hydro Storage')
            %don't plot state-of charge on electric plot
        else
            stor(end+1) = i;
        end
        dataI = Data(:,i);
        include = true;
    elseif isfield(Plant.Generator(i).QPform.output,S) %utilities
        include = true;
        dataI = Data(:,i)*Plant.Generator(i).QPform.output.(S)(1); %Hratio for CHP generators      
    end
    if include
        name(end+1,1) = {Plant.Generator(i).Name};
        Data2(:,end+1) = dataI;
        if any(stor==i)% energy storage
            storNames(end+1) ={Plant.Generator(i).Name};
        end
    end
end
end%Ends function sortForPlot

function dataPlot(h,nPlot,LINE,Data2,hours,dt,s1,name,Names,colorsPlot,tSize,S,dateText)
nnList = {'Electrical';'DistrictHeat';'DistrictCool';'Hydro';'DirectCurrent';'CoolingWater';'Transmission1';'Hydrogen';'LiqHydrogen';'Heating2';'BuildingTemp'};
LabelYaxis = {'Electrical Generation ';'Heating Generation ';'Cooling Generation ';'Withdrawls ';'DC Electrical Generation ';'CoolingWater Temperature ';'230kV Transmission';'Hydrogen ';'Liquid H2 ';'Steam ';'Temperature'};
index = nonzeros((1:11)'.*strcmp(S,nnList));

if tSize==12
    s0 = 1;
else
    s0 = s1;
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

if ~LINE && any(any(Data3<0))
    pos = sum(max(Data2(s0:end-1,:),0),2);
    neg = sum(min(0,Data2(s0:end-1,:)),2);
    OoM = max(1,.2+log10(max(pos)-min(neg)));
else
    OoM = max(1,.2+log10(max(sum(Data2(s0:end,:),2))));
end
ylab = LabelYaxis{index};
if OoM>6.30103
    if index == 4
        units = '1e9 cfs';
    else
        units = 'GW';
    end
    OoM = OoM-6;
    Data3 = Data3/1e6;
elseif OoM>3.30103
    if index == 4
        units = '1e6 cfs';
    else
        units = 'MW';
    end
    OoM = OoM-3;
    Data3 = Data3/1000;
else
    if index == 4
        units = '1000 cfs';
    elseif index == 6 || index == 11
        units = 'F';
    else
        units = 'kW';
    end
end
if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
    Yspace = 10^(OoM-1);
    Ymax = 10^OoM;
elseif (OoM-floor(OoM))> 0.6990 
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

if LINE
    str2 = 'Color';
    h1 = plot(h,plotTime,Data3,'LineWidth',3);
else
    posBars = max(0,Data3);
    negBars = min(0,Data3);
    Ikeep = any(negBars<0);
    negBars = negBars(:,Ikeep);
    negBarNames = name(Ikeep);

    str2 = 'FaceColor';
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

if tSize==12 && s1>1
    plot(h,[hours(s1),hours(s1)],[Ymin,Ymax],'c--')   
end
xlim(h,[hours(s0), hours(end)])
ylim(h,[Ymin,Ymax])
axTickY = Ymin:Yspace:Ymax;
set(h,'YTick',axTickY,'YTickLabel', {axTickY},'FontSize',tSize-2)
set(h,'XTick',axTick,'XTickLabel', {axIndex})
if nPlot==1
    xlabel(h,dateText{1},'Color','k','FontSize',tSize)
else
    xlabel(h,dateText{2},'Color','k','FontSize',tSize)
end
ylabel(h,strcat(ylab,{'  ('},units,')'),'Color','k','FontSize',tSize)
end%Ends function dataPlot

function PlotWithStorage(Data,StorageState,dt,hours,s1,colorsPlot,storNames,Names,tSize,stor,name,nPlot,handles,S,dateText)
if tSize==12
    s0 = 1;
else
    s0 = s1;
end
hName = handles.(strcat('ResultName',num2str(nPlot)));
h = handles.(strcat('ResultPlot',num2str(nPlot)));
h2 = handles.(strcat('ResultPlot',num2str(nPlot),'b'));
if length(stor)>2
    if strcmp(get(handles.(strcat('ResultStorageSelect',num2str(nPlot))),'Visible'),'off')
        set(handles.(strcat('ResultStorageSelect',num2str(nPlot))),'Visible','on','value',1,'string',storNames)
    end
    I = true(1,length(Data(1,:)));
    I(stor) = false;
    storNum = get(handles.(strcat('ResultStorageSelect',num2str(nPlot))),'value');%select only 1 storage at a time
    stor = stor(storNum);
    I(stor) = true;
    storNames = storNames(storNum);
    name = name(I);
    Data = Data(:,I);
    StorageState = StorageState(:,I);
else
    if strcmp(get(handles.(strcat('ResultStorageSelect',num2str(nPlot))),'Visible'),'on')
        set(handles.(strcat('ResultStorageSelect',num2str(nPlot))),'Visible','off')
    end
end
dataPlot(h,nPlot,false,Data,hours,dt,s1,name,Names,colorsPlot,tSize,S,dateText)
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
end%Ends function PlotWithStorage