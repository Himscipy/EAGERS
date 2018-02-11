function PlotHistorical(handles,demand,tab)
global TestData
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end

if strcmp(get(handles.sliderZoom1,'Visible'),'off')
    set(handles.(strcat('sliderZoom',num2str(tab))),'Visible','on','value',1);set(handles.(strcat('sliderDate',num2str(tab))),'Visible','on','value',0);
    set(handles.(strcat('textDay',num2str(tab))),'Visible','on'); set(handles.(strcat('textAllData',num2str(tab))),'Visible','on'); set(handles.(strcat('textHorizon',num2str(tab))),'Visible','on');
end
%find the current date
DateSim = TestData.Timestamp(1) + get(handles.(strcat('sliderDate',num2str(tab))),'Value');
maxSlider = get(handles.(strcat('sliderZoom',num2str(tab))),'Max');
if get(handles.(strcat('sliderZoom',num2str(tab))),'Value')==maxSlider
    DateEnd = TestData.Timestamp(end);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<1
    DateEnd = min(TestData.Timestamp(end),DateSim + 1);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<2
    DateEnd = min(TestData.Timestamp(end),DateSim + 7);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<3
    DateEnd = min(TestData.Timestamp(end),DateSim + 31);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<4
    DateEnd = min(TestData.Timestamp(end),DateSim + 365);
end
Xi = nnz(TestData.Timestamp<=DateSim);
Xf = nnz(TestData.Timestamp<=DateEnd);
Date = TestData.Timestamp(Xi:Xf);
if strcmp(demand,'Electrical Demand')
    data = TestData.Demand.E(Xi:Xf,:);
elseif strcmp(demand,'Cooling Demand')
    data = TestData.Demand.C(Xi:Xf,:);
elseif strcmp(demand,'Heating Demand')
    data = TestData.Demand.H(Xi:Xf,:);
elseif strcmp(demand,'Water Demand')
    data = TestData.Demand.W(Xi:Xf,:);
else
    data = TestData.Building.(demand)(Xi:Xf,:);
end

units = ' (kW)';
color = {'black'};
if strcmp(demand,'Cooling Demand')
    color = {'blue'};
elseif strcmp(demand,'Heating Demand')
    color = {'red'};
elseif strcmp(demand,'Water Demand')
    color = {'cyan'};
    units = ' (1000 CFS)';
end
if tab ==1 
    h = handles.axesMain;
    PlotData(h,Date,data,[demand, units],color)
    h = handles.axesCumulative;
    PlotHistogram(h,data,[demand, units])
elseif tab ==2
    h = handles.axesBuildLoad;
    PlotData(h,Date,data,[demand, units],color)
    h = handles.axesBuildCumulative;
    PlotHistogram(h,data,[demand, units])
end
end%Ends function PlotHistorical

function PlotData(h,TimeVec,Y,Ylab,color)
cla(h);
D = datevec(TimeVec(1));
days = max(1,round(TimeVec(end)-TimeVec(1)));
dNum = TimeVec(1);
monthvec = [0 31 59 90 120 151 181 212 243 273 304 334 365];
leapyear = mod((D(1)-2004),4)==0;%if it is a leap year it will have a remainder of 0
if leapyear
    monthvec(3:end) = monthvec(3:end)+1;
end
if days>sum(monthvec) %if you are plotting multiple years make the ticks in years
    E = datevec(TimeVec(end));
    plotAxis = dNum+linspace(0,round(TimeVec(end)-TimeVec(1)),E(1)-D(1)+1);
elseif days==sum(monthvec) %if you are plotting a year, make the ticks in months
    plotAxis = dNum+[monthvec(D(2):end),monthvec(12)+monthvec(1:D(2))]-monthvec(D(2))-D(3);%make sure to add the right number of days given the month and day of the year.
elseif days > 31
    [y,m1,~] = datevec(TimeVec(1));
    [~,m2,~] = datevec(TimeVec(end));
    M=0;
    for i=1:1:(m2-m1)
        d = datenum([y,m1+i,1])-datenum([y,m1+i-1,1]);%days in month
        M(end+1) = M(end) + d;
    end
    plotAxis = dNum+M;
elseif days>1
    plotAxis = dNum+linspace(0,days,days+1);
else
    hours = floor(24*(TimeVec(end)-TimeVec(1)));
    plotAxis = dNum+linspace(0,hours/24,floor(hours/3)+1);  
end

[m,n] = size(Y);
for i = 1:1:n
    plot(h,TimeVec,Y(:,i),color{i});
end
ylabel(h,Ylab)
xlim(h,[TimeVec(1) TimeVec(end)])
set(h,'xtick',plotAxis)
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if days>366
    datetick(h,'x','yyyy','keeplimits')
    xlabel(h,'Year') 
elseif days==365|| days==366
    datetick(h,'x','mmmdd','keeplimits')
    xlabel(h,num2str(D(1)))
elseif days>=28 && days<=31
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days==7
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days ==1
    datetick(h,'x','HH','keeplimits','keepticks')
    xlabel(h,strcat(['Hours of ', monthLabel(D(2),:), num2str(D(3))]))
end
end %Ends function PlotData

function PlotHistogram(h,Data,Xlab)
cla(h);
N = length(Data);
dSort = sort(Data);
Xi = max(1,floor(0.01*N));
Xf = ceil(0.99*N);
range = dSort(Xf)-dSort(Xi);
OoM = log10(range);
if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^(OoM-1);
elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^floor(OoM);
elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
    Xspace = .5*10^floor(OoM);
else  %count in increments of 2, 20, 200 or 2000 etc
    Xspace = .2*10^floor(OoM);
end
hist = [];
label = {};
Xi = ceil(dSort(Xi)/Xspace)*Xspace;
hist(end+1) = nnz(dSort<=Xi);
label(end+1) = {strcat('<',num2str(Xi))};
while Xi<dSort(Xf)
    Xi = Xi + Xspace;
    hist(end+1) = (nnz(dSort<=Xi) - sum(hist(1:end)));
    label(end+1) = {strcat(num2str(Xi-Xspace),'--',num2str(Xi))};
end
hist = hist/N*100;
bar(h,hist)
ylabel(h,'Percent of Time Within Range')
% set(h,'XTickLabel', label)
set(h,'XTickLabel','')
xlim(h,[0.5,length(hist)+.5]);
pos = get(h,'position');
t = text(0,0,Xlab);
set(t,'Parent',h,'Units','characters','HorizontalAlignment','center','Position',[pos(3)/2,-6,0]);
% Place the text labels
for i = 1:length(hist)
    t = text(0,0,label{i});
    xpos = i*pos(3)/length(hist)- 0.5*pos(3)/length(hist);
    set(t,'Parent',h,'Units','characters','HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'Position',[xpos,0,0]);
end
end %Ends function PlotHistogram