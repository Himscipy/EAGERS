function varargout = DISPATCH(varargin)
% DISPATCH MATLAB code for DISPATCH.fig
%      DISPATCH, by itself, creates a new DISPATCH or raises the existing
%      singleton*.
%
%      H = DISPATCH returns the handle to a new DISPATCH or the handle to
%      the existing singleton*.
%
%      DISPATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPATCH.M with the given input arguments.
%
%      DISPATCH('Property','Value',...) creates a new DISPATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DISPATCH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DISPATCH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
 
% Edit the above text to modify the response to help DISPATCH
 
% Last Modified by GUIDE v2.5 24-May-2017 20:29:02
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DISPATCH_OpeningFcn, ...
                   'gui_OutputFcn',  @DISPATCH_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
 
% --- Executes just before DISPATCH is made visible.
function DISPATCH_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for EAGERS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant GENINDEX TestData RealTimeData
GENINDEX = 1;
set(gcf,'Name','DISPATCH')
movegui(gcf,'center');
if strcmp(Plant.optimoptions.method,'Planning')
    Plant.optimoptions.method = 'Dispatch';
end
TestData = Plant.Data;%% Revise this so you can pull from more than what is loaded in Plant
RealTimeData = Plant.Data; %this will get replaced later, but needed for forecasting initially
%% Set up Main Tabs
% Assumptions:
% 1. Tags of main tab static text boxes are of form, 'MainTab1',
% 'MainTab2', etc.
% 2. Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2',
% etc.
TabText = {'Main Window';'Market Services';'Historian/Forecast';'Network View';'Settings'};
% set(hObject,'UserData',TabText);
for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'Units','characters','String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'Units','characters','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        pos = get(handles.MainTab1','Position');
        set(handles.MainTab1,'Position',[pos(1),pos(2),pos(3),pos(4)+.5])
    else
        pos2 = get(handles.(strcat('uipanelMain',j)),'Position');
        set(handles.(strcat('uipanelMain',j)),'Position',[pan1pos(1), pan1pos(2)+(pan1pos(4)-pos2(4)), pos2(3), pos2(4)])
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end

%make gen list for tabs 1 and 5
list={};
for i=1:length(Plant.Generator)
    if Plant.Generator(i).Enabled ==1
        list(end+1) = {Plant.Generator(i).Name};
    elseif Plant.Generator(i).Enabled ==0
        list(end+1) = {Plant.Generator(i).Name};
    end 
end
set(handles.uipanelMain1,'UserData',list)

networkNames = fieldnames(Plant.Network);
handles = PlotAxes_Setup(hObject, eventdata, handles,networkNames,1);
handles = PlotAxes_Setup(hObject, eventdata, handles,networkNames,3);

set(handles.popupmenuColor,'String',{'parula';'autumn';'cool';'spring';'summer';'winter';},'Value',1);

colormap(handles.ResultPlot1,'parula');
set(hObject,'UserData',colormap(handles.ResultPlot1));

handles = GenList_Make(handles,1);
handles = GenList_Make(handles,3);
handles = GenList_Make(handles,5);

%Set handles for Settingss
set(handles.editCommandOnOff,'string','--')
set(handles.editCommandSet,'string','--')
set(handles.editMeasureOnOff,'string','--')
set(handles.editMeasureInput,'string','--')
set(handles.editMeasurePrimary,'string','--')
set(handles.editMeasureSecondary,'string','--')

set(handles.constant, 'value', strcmp(Plant.optimoptions.tspacing,'constant'));
set(handles.linear, 'value', strcmp(Plant.optimoptions.tspacing, 'linear'));
set(handles.logarithm, 'value', strcmp(Plant.optimoptions.tspacing, 'logarithm'));
set(handles.manual, 'value', strcmp(Plant.optimoptions.tspacing, 'manual'));
set(handles.Interval,'string',Plant.optimoptions.Interval);
set(handles.Horizon, 'string', Plant.optimoptions.Horizon);
set(handles.Resolution, 'string', Plant.optimoptions.Resolution);
set(handles.scaletime, 'string', Plant.optimoptions.scaletime);

set(handles.sequential, 'value', Plant.optimoptions.sequential);
set(handles.simultaneous, 'value', ~Plant.optimoptions.sequential);
set(handles.excessHeat, 'value', Plant.optimoptions.excessHeat);
set(handles.nsSmooth, 'string', Plant.optimoptions.nsSmooth);

set(handles.NoMixedInteger, 'value', ~Plant.optimoptions.MixedInteger);
set(handles.MixedInteger, 'value', Plant.optimoptions.MixedInteger);
set(handles.NREL, 'value', 0);

set(handles.noSpinReserve, 'value', ~Plant.optimoptions.SpinReserve);
set(handles.SpinReserve, 'value', Plant.optimoptions.SpinReserve);
set(handles.SpinReservePerc, 'string', Plant.optimoptions.SpinReservePerc);
set(handles.editBuffer, 'string', Plant.optimoptions.Buffer);

set(handles.Control, 'value', strcmp(Plant.optimoptions.method,'Control'));
set(handles.Dispatch, 'value', strcmp(Plant.optimoptions.method,'Dispatch'));
set(handles.Topt, 'string', Plant.optimoptions.Topt);
set(handles.Tmpc, 'string', Plant.optimoptions.Tmpc);

set(handles.SES, 'value', false);
set(handles.ARIMA, 'value', false);
set(handles.NeuralNet, 'value', false);
set(handles.Surface, 'value', false);
set(handles.Perfect, 'value', false);
set(handles.(Plant.optimoptions.forecast), 'value', true);

%Update forecastmethod handles
days = ceil(Plant.Data.Timestamp(end) - Plant.Data.Timestamp(1));
set(handles.sliderZoom,'Min',0,'Max',1,'Value',1/7)
set(handles.sliderDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/(days-1),10/(days-1)])
set(handles.sliderStartDate,'Min',0,'Max',days,'Value',0,'SliderStep',[1/(days-1),10/(days-1)])
%%put something into axes
sliderStartDate_Callback(hObject, eventdata, handles)
LineGraph_Callback(hObject, eventdata, handles)

% --- Outputs from this function are returned to the command line.
function varargout = DISPATCH_OutputFcn(hObject, eventdata, handles) 
%ask user to save?

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
global GENINDEX
n = get(hObject,'Tag');
n = n(end);
m = [];
i = 1;
% Find out which tab is currently selected
while isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        m = i;
    else i = i+1;
    end
end
m = num2str(m);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% Change color
bColor = get(handles.(strcat('MainTab',m)),'BackgroundColor');
set(handles.(strcat('MainTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% Change dimensions
pos = get(handles.(strcat('MainTab',m)),'Position');
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.5])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.5])

% Change visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')

%% actions specific to one tab or another
if strcmp(n,'3')
    GENINDEX = [];
    set(handles.Demands,'Visible','off')
    ForecastPlot
elseif strcmp(n,'1')
    sliderStartDate_Callback(hObject, eventdata, handles)
elseif strcmp(n,'4')

end
 
% --- Executes on button press in Switch.
function Switch_Callback(hObject, eventdata, handles)
%%send user back to EPT, pass along the plant generators 
close
MainScreen1

% --- Creates plots on the Main Window and Forecasting Tabs
function handles = PlotAxes_Setup(hObject, eventdata, handles,networkNames,tab)
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
nPlot = length(networkNames);
if nPlot ==1 %1 very large plot
    Pos = [52 6 150 30];
else %1 large plot and n-1 smaller plots
    Pos = [52 10 80 26;];
    for j = 2:nPlot
        Pos(j,:) = [153 3+(j-2)*40/(nPlot-1) 50 min(20,(40/(nPlot-1)-6))];
    end
end
if tab==1
    name = 'Result';
else name = 'Forecast';
end
quote='''';
for i = 1:1:nPlot
    %primary axes
    handles.(strcat(name,'Plot',num2str(i))) = axes('Units','characters',...
        'Position', Pos(i,:),'NextPlot','add',...
        'Tag', strcat(name,'Plot',num2str(i)),...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'Visible','on');
    %secondary axes (y-label on right)
    handles.(strcat(name,'Plot',num2str(i),'b')) = axes('Units','characters',...
        'Position', Pos(i,:),'NextPlot','add','color','none',...
        'Tag', strcat(name,'Plot',num2str(i),'b'),...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'YAxisLocation','right',...
        'Visible','on');
    callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'PlotSwitch',quote,',hObject)');
    handles.(strcat(name,'Name',num2str(i))) = uicontrol('Style', 'pushbutton', 'String', networkNames{i},...
        'Units','characters',...
        'Position', [Pos(i,1)+Pos(i,3)/10,Pos(i,2)+Pos(i,4)+.25,35,2],...
        'BackGroundColor',[0.94 0.94 0.94],...
        'Tag', strcat(name,'Name',num2str(i)),...
        'FontSize',16,...
        'FontUnits','points',...
        'FontWeight','bold',...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'Callback',eval(callback),...
        'Visible','on',...
        'UserData',i);
    if tab==1
        handles.(strcat(name,'StorageSelect',num2str(i))) = uicontrol('Style', 'popupmenu', 'String', networkNames{i},...
            'Units','characters',...
            'Position', [Pos(i,1)+40,Pos(i,2)+Pos(i,4)+.25,30,2],...
            'Tag', strcat(name,'StorageSelect',num2str(i)),...
            'FontSize',11,...
            'FontUnits','points',...
            'FontWeight','normal',...
            'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
            'Visible','off');
    end
end
%% add pop-up menue with storage names if more than 2 storage devices of this type.

function PlotSwitch(hObject)
%if you pick the main plot and it is network1, swap with plot 2, 
% otherwise swap the small plot picked and main plot
global Plant
handles = guihandles;
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
nPlot = length(networkNames);
net = get(hObject,'String');
plotI = get(hObject,'UserData');
netI = find(strcmp(net,networkNames));
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    screen = 'Result';
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    screen = 'Forecast';
end
for i = 1:1:nPlot %put names back to default pos
    set(handles.(strcat(screen,'Name',num2str(i))),'String',networkNames{i});
end
    
if netI == plotI %graph is currently in its default position (move to primary)
    if netI == 1 && nPlot>1%already is primary swap with 2
        set(handles.(strcat(screen,'Name1')),'String',networkNames{2});
        set(handles.(strcat(screen,'Name2')),'String',networkNames{1});
    else
        set(handles.(strcat(screen,'Name1')),'String',networkNames{netI});
        set(handles.(strcat(screen,'Name',num2str(netI))),'String',networkNames{1});
    end
end
if strcmp(screen,'Forecast')
    ForecastPlot
elseif strcmp(screen,'Result')
    A = get(handles.ResultPlot1,'UserData');
    ForecastTime = A{1};
    Forecast = A{2};
    HistoryTime = A{3};
    History = A{4};
    plotDispatch(handles,ForecastTime,Forecast,HistoryTime,History)
end

function handles = GenList_Make(handles,tab)
global Plant Model_dir
quote='''';
list = get(handles.uipanelMain1,'UserData');
if tab == 1
    r = 'Main';
elseif tab == 3
    r = 'Plot';
elseif tab == 5
    r = 'Comm';
end
if length(list)>12
    set(handles.(strcat('NextGen',num2str(tab))),'Visible','on','UserData',1)
    set(handles.(strcat('PrevGen',num2str(tab))),'UserData',1)
end
%Makes buttons for GUI && status buttons next to corresponding generator
nG = length(list);
colorVec = get(gcf,'UserData');
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG));
for i=1:1:nG 
    num = num2str(i);
    curtab = floor((i-1)/12)+1;
    prev = 12*(curtab-1);
    if curtab==1
        vis = 'on';
    else vis = 'off';
    end
    pos = 45 - 2*(i-prev);
    callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'Gen',r,'_Callback',quote,',hObject,eventdata,guidata(hObject))');
    handles.(strcat('Generator',num2str(tab),'_',num)) = uicontrol('Style', 'pushbutton', 'String', list{i},...
    'Units','characters',...
    'Position', [1 pos 25 1.8],...
    'Tag', strcat('Generator',num2str(tab),'_',num),...
    'FontSize', 10,...
    'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
    'Callback',eval(callback),...
    'Visible',vis,...
    'UserData',i,...
    'BackgroundColor',colorsPlot(i,:));
    if tab ==1 %Only make Status buttons on Main Window
        pos = 45.5 - 2*(i-prev);
        callback = strcat('@(hObject,eventdata)DISPATCH(',quote,'Status_Callback',quote,',hObject,eventdata,guidata(hObject))');
        if Plant.Generator(i).Enabled
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
        else
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
        end
        if license('test','Image Processing Toolbox')
            set(handles.Switch,'Units','pixels');
            pos1 = get(handles.Switch,'Position');
            set(handles.Switch,'Units','characters');
            pos2 = get(handles.Switch,'Position');
            x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
        end
        if Plant.Generator(i).Enabled
            enableGen  = 'bold';
        else enableGen  = 'normal';
        end
        handles.(strcat('GeneratorStat',num2str(tab),'_',num)) = uicontrol('Style', 'pushbutton', 'String', '',...
        'Units','characters',...
        'Position', [27 pos 3 1],...
        'Tag', strcat('GeneratorStat',num2str(tab),'_',num),...
        'cdata', x,...
        'FontWeight',enableGen,...
        'Parent', handles.uipanelMain1,...
        'Callback',eval(callback),...
        'Visible',vis,...
        'UserData',i);
    end
end
% --- Executes on button press in PrevGen1.
function PrevGen_Callback(hObject, eventdata, handles)
handles = guihandles;
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    tab=1;
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    tab=3;
elseif strcmp(get(handles.uipanelMain5,'Visible'),'on')
    tab=5;
end
list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.(strcat('PrevGen',num2str(tab))),'UserData');%current page of the list
if page<2
    set(handles.(strcat('PrevGen',num2str(tab))),'Visible','off','UserData',page-1)
else
    set(handles.(strcat('PrevGen',num2str(tab))),'Visible','on','UserData',page-1);
end
set(handles.(strcat('NextGen',num2str(tab))),'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('Generator',num2str(tab),'_',j)),'Visible','off');
         if tab==1
             set(handles.(strcat('GeneratorStat',num2str(tab),'_',j)),'Visible','off');
         end
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('Generator',num2str(tab),'_',j)),'Visible','on');
     if tab==1
         set(handles.(strcat('GeneratorStat',num2str(tab),'_',j)),'Visible','on');
     end
end


% --- Executes on button press in NextGen1.
function NextGen_Callback(hObject, eventdata, handles)
handles = guihandles;
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    tab=1;
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    tab=3;
elseif strcmp(get(handles.uipanelMain5,'Visible'),'on')
    tab=5;
end

list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.(strcat('PrevGen',num2str(tab))),'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.(strcat('NextGen',num2str(tab))),'Visible','off','UserData',page+1)
else
    set(handles.(strcat('NextGen',num2str(tab))),'Visible','on','UserData',page+1);
end
set(handles.(strcat('PrevGen',num2str(tab))),'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('Generator',num2str(tab),'_',j)),'Visible','off');
     if tab==1
         set(handles.(strcat('GeneratorStat',num2str(tab),'_',j)),'Visible','off');
     end
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('Generator',num2str(tab),'_',j)),'Visible','on');
         if tab==1
             set(handles.(strcat('GeneratorStat',num2str(tab),'_',j)),'Visible','on');
         end
    end
end

%When Component Buttons on the main tab are clicked
function GenMain_Callback(hObject, eventdata, handles)
global Plant GENINDEX
gen = get(hObject,'String');
GENINDEX = get(hObject,'UserData');
size = num2str(Plant.Generator(GENINDEX).Size);
set(handles.SelGen,'Title',Plant.Generator(GENINDEX).Name,'UserData',GENINDEX)
if ~isempty(strfind(gen,'Utility'))
    set(handles.GenSpec1,'String','Inf')
else
    set(handles.GenSpec1,'String',size)
end
if Plant.Generator(GENINDEX).Enabled == 1
    set(handles.GenEnable,'Value',1)
    set(handles.GenDisable,'Value',0)
elseif Plant.Generator(GENINDEX).Enabled == 0
    set(handles.GenEnable,'Value',0)
    set(handles.GenDisable,'Value',1)
end
 

%When Status colors are clicked
function Status_Callback(hObject, eventdata,handles)
global Plant Model_dir
i = get(hObject,'UserData');
size = num2str(Plant.Generator(i).Size);
if ~isempty(strfind(Plant.Generator(i).Name,'Utility'))
    set(handles.GenSpec1,'String','Inf')
else
    set(handles.GenSpec1,'String',size)
end
set(handles.SelGen,'Title',Plant.Generator(i).Name,'UserData',i)
if strcmp(get(hObject,'FontWeight'),'bold')
    Plant.Generator(i).Enabled = 0;
    set(hObject,'FontWeight','normal')
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
else
    Plant.Generator(i).Enabled = 1;
    set(hObject,'FontWeight','bold')
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
end
set(handles.GenEnable,'value',Plant.Generator(i).Enabled)
set(handles.GenDisable,'Value',~Plant.Generator(i).Enabled)
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
set(hObject,'cdata',x)

% --- Executes on slider movement.
function sliderStartDate_Callback(hObject, eventdata, handles)
global Plant Last24hour DateSim 
handles = guihandles;
%find the current date
DateSim = Plant.Data.Timestamp(1) + get(handles.sliderStartDate,'Value');
Last24hour = [];
if Plant.Data.Timestamp(1)>(DateSim-1) && Plant.Data.Timestamp(end)<DateSim
    Timestamp = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
else %need to have this in terms of the first timestep
    Timestamp = DateSim+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
end
Last24hour = GetHistoricalData(Timestamp);
Last24hour.Timestamp = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));

Data = GetHistoricalData(DateSim); 
if isfield(Plant,'Dispatch') && ~isempty(Plant.Dispatch) && isfield(Plant.Dispatch,'Timestamp') && ~isempty(Plant.Dispatch.Timestamp)
    I = nnz(Plant.Dispatch.Timestamp<=DateSim & Plant.Dispatch.Timestamp>0);
else I = 0;
end
if I == 0 || Plant.Dispatch.Timestamp(I)<(DateSim - Plant.optimoptions.Resolution/24)%if it has not been run at this timestep
    if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
        initializeOptimization
    end
    automaticInitialCondition(Data);
    ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
    Forecast = updateForecast(ForecastTime(2:end),Data);
    Dispatch = DispatchLoop(ForecastTime,Forecast,[]);
    for i = 1:1:length(Plant.Generator)
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Dispatch(1,i) = RenewableOutput(Plant.Generator(i).VariableStruct,DateSim,'Actual');
        end
    end
    Forecast = Dispatch;
    History = [];
    HistoryTime = [];
else
    if Plant.Dispatch.Timestamp(I+1)>0 && (Plant.Dispatch.Timestamp(I+1)-DateSim)<(DateSim-Plant.Dispatch.Timestamp(I))
        I = I+1; %The next time step is actually closer
    end
    ForecastTime = Plant.Predicted.Timestamp(2:end,I);
    Forecast = Plant.Predicted.GenDisp(2:end,:,I);
    backSteps = min(I,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
    History = Plant.Dispatch.GeneratorState(I-backSteps+1:I,:);
    HistoryTime = Plant.Dispatch.Timestamp(I-backSteps+1:I);
end
plotDispatch(handles,ForecastTime,Forecast,HistoryTime,History)

function sliderStartDate_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
global Plant Virtual RealTime DispatchWaitbar DateSim Last24hour GenAvailTime RestartTime RealTimeData
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
if get(handles.VirtualMode,'Value') ==1 
    Virtual = 1;
    RealTime = 0;
    Plant.optimoptions.fastsimulation = 1;
elseif get(handles.ObserverMode,'Value') == 1
    Virtual = 1;
    RealTime = 0;
    Plant.optimoptions.fastsimulation = 0;
elseif get(handles.ControllerMode,'Value') == 1
    Virtual = 0;
    RealTime = 1;
    Plant.optimoptions.fastsimulation = 0;
end
set(handles.Start,'Value',1);%reset start button
set(handles.Stop,'Value',0);%reset stop button

%%Need to select a starting Date
DateSim = Plant.Data.Timestamp(1) + get(handles.sliderStartDate,'Value');
if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
    initializeOptimization%Load generators, build QP matrices
end
% clear & initialize variables
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<=DateSim+Plant.optimoptions.Interval+Plant.optimoptions.Horizon/24);
if any(strcmp(Plant.optimoptions.method,{'Dispatch';'Planning'}))
    RealTimeData = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);%create test data at correct frequency
else %actual control
    nG = length(Plant.Generator);
    for i = 1:1:nG
        if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
            RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
        else
            RestartTime(i) = 0;
        end
    end
    GenAvailTime = ones(1,nG).*DateSim; %  Global vairables needed in controller mode
    RealTimeData = interpolateData(Plant.optimoptions.Tmpc,Xi,Xf,0.00);%create test data at correct frequency
end

nG = length(Plant.Generator);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG;
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution;
Plant.Dispatch.Temperature = zeros(NumSteps,1);
Plant.Dispatch.Timestamp = zeros(NumSteps,1);
Plant.Dispatch.GeneratorState = zeros(NumSteps,nG+nL);
Plant.Predicted.GenDisp = zeros(round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution)+1,nG+nL,NumSteps);
Plant.Predicted.Timestamp = zeros(round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution)+1,NumSteps);
Outs =  fieldnames(Plant.Data.Demand);
for i = 1:1:length(Outs)
    loads = length(Plant.Data.Demand.(Outs{i})(1,:));
    Plant.Dispatch.Demand.(Outs{i}) = zeros(NumSteps,loads);
    Plant.RunData.Demand.(Outs{i}) = zeros(NumSteps,loads);
end

Last24hour = [];
TimeYesterday = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
if Plant.Data.Timestamp(1)<(DateSim-1) && Plant.Data.Timestamp(end)>=DateSim
    Last24hour = GetCurrentData(TimeYesterday);
else %need to have this in terms of the first timestep
    Last24hour = GetCurrentData(TimeYesterday+1);
    Last24hour.Timestamp = TimeYesterday;
end

% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;
if K ==1
    Plant.Dispatch.GeneratorState(1,:) = manualInitialCondition;
else
    Plant.Dispatch.GeneratorState(1,:) = automaticInitialCondition(GetCurrentData(DateSim));
end
Plant.Dispatch.Timestamp(1) = DateSim;
Si=1; %counter for # of times dispatch loop has run
LastDispatch = [];
mainFig = gcf;
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
timers = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
while Si<NumSteps
    Date = DateSim+[0;Time/24];
    Data = GetCurrentData(DateSim);
    Forecast = updateForecast(Date(2:end),Data);%% function that creates demand vector with time intervals coresponding to those selected
    [LastDispatch,timers(Si,:)] = DispatchLoop(Date,Forecast,LastDispatch);
%     disp(strcat('FistDisp:',num2str(timers(Si,1))));
%     disp(strcat('StebByStep:',num2str(timers(Si,2))));
%     disp(strcat('FinalDisp:',num2str(timers(Si,3))));
    for i = 1:1:length(Plant.Generator)
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            LastDispatch(1,i) = RenewableOutput(Plant.Generator(i).VariableStruct,DateSim,'Actual');
        end
    end
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end
    backSteps = min(Si,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
    History = Plant.Dispatch.GeneratorState(Si-backSteps+1:Si,:);
    HistoryTime = Plant.Dispatch.Timestamp(Si-backSteps+1:Si,:);
    
    handles = guihandles(mainFig);
    updateGUIstatus(handles,LastDispatch(1:2,:),History)
    plotDispatch(handles,Date(2:end),LastDispatch(2:end,:),HistoryTime,History)
    
    Si = StepDispatchForward(Si,Date,Data,Forecast,LastDispatch);
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
end
if strcmp(Plant.optimoptions.method,'Dispatch')
    Plant.Cost.Dispatch = NetCostCalc(Plant.Dispatch.GeneratorState,Plant.Dispatch.Timestamp,'Dispatch');
elseif strcmp(Plant.optimoptions.method,'Control')
    Plant.Cost.RunData = NetCostCalc(Plant.RunData.GeneratorInput,Plant.RunData.Timestamp,'Input');
end
% Plant.Baseline = RunBaseline(Plant.Dispatch.GeneratorState(1,:)); %finish simulation by running baseline
Stop_Callback(hObject, eventdata, handles)

% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
%stops dispatch
global DispatchWaitbar Virtual
set(handles.Stop,'Value',1);%reset stop button
set(handles.Start,'Value',0);%reset start button
Virtual = 0;
close(DispatchWaitbar)
DispatchWaitbar=[];

% --- Executes on button press in GenEnable.
function GenEnable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 1;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',x)

% --- Executes on button press in GenDisable.
function GenDisable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 0;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
if license('test','Image Processing Toolbox')
    set(handles.Switch,'Units','pixels');
    pos1 = get(handles.Switch,'Position');
    set(handles.Switch,'Units','characters');
    pos2 = get(handles.Switch,'Position');
    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
end
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',x)
        
%% for manual control?
function GenStatus1_Callback(hObject, eventdata, handles)
function GenStatus1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenStatus2_Callback(hObject, eventdata, handles)
function GenStatus2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec1_Callback(hObject, eventdata, handles)
function GenSpec1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec2_Callback(hObject, eventdata, handles)
function GenSpec2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LineGraph.
function LineGraph_Callback(hObject, eventdata, handles)
if get(handles.StackedGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.LineGraph,'Value',1); %was already pressed
end
handles = guihandles;
A = get(handles.ResultPlot1,'UserData');
ForecastTime = A{1};
Forecast = A{2};
HistoryTime = A{3};
History = A{4};
plotDispatch(handles,ForecastTime,Forecast,HistoryTime,History)

% --- Executes on button press in StackedGraph.
function StackedGraph_Callback(hObject, eventdata, handles)
if get(handles.LineGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.StackedGraph,'Value',1); %was already pressed
end
handles = guihandles;
A = get(handles.ResultPlot1,'UserData');
ForecastTime = A{1};
Forecast = A{2};
HistoryTime = A{3};
History = A{4};
plotDispatch(handles,ForecastTime,Forecast,HistoryTime,History)

% --- Executes on button press in AutoControl.
function AutoControl_Callback(hObject, eventdata, handles)
if get(handles.AutoControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.AutoControl,'Value',1); %was already pressed
end

% --- Executes on button press in ManualControl.
function ManualControl_Callback(hObject, eventdata, handles)
if get(handles.ManualControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ManualControl,'Value',1); %was already pressed
end

% --- Executes on button press in VirtualMode.
function VirtualMode_Callback(hObject, eventdata, handles)
if get(handles.VirtualMode,'Value')==1
    a = get(handles.VirtualMode,'BackgroundColor');
    c = get(handles.VirtualMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.VirtualMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.VirtualMode,'Value',1); %was already pressed
end

% --- Executes on button press in ObserverMode.
function ObserverMode_Callback(hObject, eventdata, handles)
if get(handles.ObserverMode,'Value')==1
    a = get(handles.ObserverMode,'BackgroundColor');
    c = get(handles.ObserverMode,'ForegroundColor');
    if get(handles.VirtualMode,'Value')==1
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ObserverMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ObserverMode,'Value',1); %was already pressed
end

% --- Executes on button press in ControllerMode.
function ControllerMode_Callback(hObject, eventdata, handles)
if get(handles.ControllerMode,'Value')==1
    a = get(handles.ControllerMode,'BackgroundColor');
    c = get(handles.ControllerMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    else
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ControllerMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ControllerMode,'Value',1); %was already pressed
end

%% Historian/ForecastMethod Tab
function GenPlot_Callback(hObject, eventdata, handles)
%plot the forecastmethod and actual dispatch of a single generator
global GENINDEX
GENINDEX = get(hObject,'UserData');
set(handles.Demands,'Visible','on')
ForecastPlot

function Demands_Callback(hObject, eventdata, handles)
%plot the forecastmethod and actual demands
global GENINDEX
GENINDEX = [];
set(handles.Demands,'Visible','off')
ForecastPlot

function ForecastPlot
global Plant Last24hour DateSim GENINDEX
handles = guihandles;
%find the current date
DateSim = Plant.Data.Timestamp(1) + get(handles.sliderDate,'Value');
Last24hour = [];
TimeYesterday = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
if Plant.Data.Timestamp(1)<(DateSim-1) && Plant.Data.Timestamp(end)>=DateSim
    Last24hour = GetHistoricalData(TimeYesterday);
else %need to have this in terms of the first timestep
    Last24hour = GetHistoricalData(TimeYesterday+1);
    Last24hour.Timestamp = TimeYesterday;
end

Data = GetHistoricalData(DateSim);        
if isempty(GENINDEX)%plot demands
    set(handles.sliderZoom,'Visible','on')
    set(handles.textHour,'Visible','on'); set(handles.textWeek,'Visible','on'); set(handles.textHorizon,'Visible','on');
    nPlot = length(fieldnames(Last24hour.Demand));
    for q = 1:1:nPlot
        set(handles.(strcat('ForecastPlot',num2str(q))),'Visible','on');
        set(handles.(strcat('ForecastPlot',num2str(q),'b')),'Visible','on');
        set(handles.(strcat('ForecastName',num2str(q))),'Visible','on');
    end
    DateEnd = DateSim + get(handles.sliderZoom,'Value')*7;
    Forecast = updateForecast(linspace(DateSim,DateEnd)',Data);%% function that creates demand vector with time intervals coresponding to those selected
    Outs =  fieldnames(Forecast.Demand);
    Xi = nnz(Plant.Data.Timestamp<(DateSim-1))+1;
    Xf = nnz(Plant.Data.Timestamp<(DateEnd) & Plant.Data.Timestamp>0);
    for i = 1:1:length(Outs)
        Actual.(Outs{i}) = Plant.Data.Demand.(Outs{i})(Xi:Xf);
    end
    ActualTime = Plant.Data.Timestamp(Xi:Xf);
    plotDemand(Forecast.Timestamp,Forecast.Demand,ActualTime,Actual)
else
    set(handles.sliderZoom,'Value',Plant.optimoptions.Horizon/24/7,'Visible','off')%hide slider
    set(handles.textHour,'Visible','off'); set(handles.textWeek,'Visible','off'); set(handles.textHorizon,'Visible','off');
    DateEnd = DateSim + get(handles.sliderZoom,'Value')*7;
    if isfield(Plant,'Dispatch')
        I = nnz(Plant.Dispatch.Timestamp<=DateSim & Plant.Dispatch.Timestamp>0);
    else I = 0;
    end
    if I == 0 || Plant.Dispatch.Timestamp(I)<(DateSim - Plant.optimoptions.Resolution/24)%if it has not been run at this timestep
        ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
        if strcmp(Plant.Generator(GENINDEX).Source,'Renewable')
            Forecast = RenewableOutput(Plant.Generator(GENINDEX).VariableStruct,ForecastTime,'Actual');
        else
            if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
                initializeOptimization
            end
            automaticInitialCondition(Data);
            Forecast = updateForecast(ForecastTime(2:end),Data);
            Dispatch = DispatchLoop(ForecastTime,Forecast,[]);
            Forecast = Dispatch(:,GENINDEX);
        end
        History = [];
        HistoryTime = [];
    else
        if Plant.Dispatch.Timestamp(I+1)>0 && (Plant.Dispatch.Timestamp(I+1)-DateSim)<(DateSim-Plant.Dispatch.Timestamp(I))
            I = I+1; %The next time step is actually closer
        end
        ForecastTime = Plant.Predicted.Timestamp(:,I);
        Forecast = Plant.Predicted.GenDisp(:,GENINDEX,I);
        backSteps = min(I,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        History = Plant.Dispatch.GeneratorState(I-backSteps+1:I,GENINDEX);
        HistoryTime = Plant.Dispatch.Timestamp(I-backSteps+1:I);
    end
    if isfield(Plant.Data,'Dispatch')
        Xi = nnz(Plant.Data.Timestamp<(DateSim-1))+1;
        Xf = nnz(Plant.Data.Timestamp<(DateEnd));
        Actual = Plant.Data.Dispatch(Xi:Xf,GENINDEX);
        ActualTime = Plant.Data.Timestamp(Xi:Xf);
    else Actual = []; ActualTime = [];
    end
    plotSchedule(ForecastTime,Forecast,HistoryTime,History,ActualTime,Actual,GENINDEX)
end

function sliderDate_Callback(hObject, eventdata, handles)
ForecastPlot
function sliderDate_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sliderZoom_Callback(hObject, eventdata, handles)
ForecastPlot
function sliderZoom_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% Network View Tab


%% Settings Tab
function removeQPmatrices
%something changed in the specifications, any previously calculated
%matrices or results do not apply
global Plant
A = {'subNet';'OpMatA';'OpMatB';'OneStep';'Online';'Dispatch';'Predicted';'RunData';'Baseline'};
for i = 1:1:length(A)
    if isfield(Plant,A{i})
        Plant.(A{i}) = [];
    end
end

function changingtimesteps_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'constant'
        Plant.optimoptions.tspacing = 'constant';
    case 'manual'
        Plant.optimoptions.tspacing = 'manual';
        prompt = {'Specify a vector of times out of one horizon for each timestep: (ex: .05,.1,.25,.75,1 would result in one timestep at 5hr, 10hr, 25hr, 75hr, and 100hr for a 100 hour horizon)'};
        dlg_title = 'Manual Timesteps';
        num_lines = 1;
        def_ans = {'0.0035,0.0070,0.0105,0.0140,0.0175,0.0210,0.0417,0.0625,0.0833,0.1042,0.125,0.1667,0.2083,0.25,0.3333,0.4167,0.5,0.5833,0.6667,0.75,0.8333,0.9167,1'};
        a = inputdlg(prompt,dlg_title,num_lines,def_ans);
        Plant.optimoptions.manualT = str2double(strsplit(a{:}, ','));
    case 'linear'
        Plant.optimoptions.tspacing = 'constant';
    case 'logarithm'
        Plant.optimoptions.tspacing = 'logarithm';
end
removeQPmatrices

function Interval_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Interval = str2double(get(handles.Interval, 'String'));
function Interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resolution_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));
removeQPmatrices

function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Horizon_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));
removeQPmatrices

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scaletime_Callback(hObject, eventdata, handles)
global Plant
%scaletime: the ratio of emulated time to time in the test-data. For example 24 hours of test data can be run in 3 hours with a scaletime of 8. scaletime enlarges any energy storage and slows down the transient response of any generators
Plant.optimoptions.scaletime = str2double(get(handles.scaletime, 'String'));
removeQPmatrices

function scaletime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChillingHeating_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'sequential'
        Plant.optimoptions.sequential = 1;
    case 'simultaneous'
        Plant.optimoptions.sequential = 0;
end
removeQPmatrices

function excessHeat_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.excessHeat = get(hObject, 'Value');
removeQPmatrices

function nsSmooth_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.nSSmooth = str2double(get(handles.nsSmooth, 'String'));
removeQPmatrices

function nsSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in simulationspeed.
function simulationspeed_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'fastsimulation'
        Plant.optimoptions.fastsimulation = 1;
    case 'slowsimulation'
        Plant.optimoptions.fastsimulation = 0;
end

% --- Executes when selected object is changed in uipanelSolverMethod.
function uipanelSolverMethod_SelectionChangedFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'NoMixedInteger'
        Plant.optimoptions.MixedInteger = 0;
    case 'MixedInteger'
        Plant.optimoptions.MixedInteger = 1;
    case 'NREL'
        %
end

% --- Executes when selected object is changed in SpinningReserve.
function SpinningReserve_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'noSpinReserve'
        Plant.optimoptions.SpinReserve = false;
    case 'SpinReserve'
        Plant.optimoptions.SpinReserve = true;
        Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
end
removeQPmatrices

function SpinReservePerc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
removeQPmatrices

function SpinReservePerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBuffer_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Buffer = str2double(get(handles.editBuffer, 'String'));
removeQPmatrices

function editBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ControlDispatch_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'Control'
        Plant.optimoptions.method = 'Control';
    case 'Dispatch'
        Plant.optimoptions.method = 'Dispatch';
end

function Topt_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Topt = str2double(get(handles.Topt, 'String'));
removeQPmatrices

function Topt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tmpc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Tmpc = str2double(get(handles.Tmpc, 'String'));
removeQPmatrices

function Tmpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NoMixedInteger_Callback(hObject, eventdata, handles)
function MixedInteger_Callback(hObject, eventdata, handles)
function noSpinReserve_Callback(hObject, eventdata, handles)
function SpinReserve_Callback(hObject, eventdata, handles)
function Control_Callback(hObject, eventdata, handles)
function Dispatch_Callback(hObject, eventdata, handles)

function ForecastMethod_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'SES'
        Plant.optimoptions.forecast = 'SNIWPE';
    case 'ARIMA'
        Plant.optimoptions.forecast = 'ARIMA';
    case 'NeuralNet'
        Plant.optimoptions.forecast = 'NeuralNet';
    case 'Surface'
        Plant.optimoptions.forecast= 'Surface';
    case 'Perfect'
        Plant.optimoptions.forecast = 'Perfect';
end

function GenComm_Callback(hObject, eventdata, handles)
global Plant
genHandles = get(handles.uipanelMain5,'Children');
j = 1;
while ~get(genHandles(j),'Value')
    j = j+1;
end
name = get(genHandles(j),'String');
i = 1;
while ~strcmp(name,Plant.Generator(i).Name)
    i = i+1;
end
set(handles.CurrentComm,'String',['Ports for: ',name])
if isfield(Plant.Generator(i).VariableStruct,'Comm')
    Comm = Plant.Generator(i).VariableStruct.Comm;
else
    Comm.OnOff = 0;
    Comm.Set = 0;
end
if isfield(Plant.Generator(i).VariableStruct,'Measure')
    Measure = Plant.Generator(i).VariableStruct.Measure;
else
    Measure.OnOff = 0;
    Measure.Input = 0;
    Measure.Electric = 0;
    Measure.Thermal = 0;
end
set(handles.editCommandOnOff,'string',num2str(Comm.OnOff))
set(handles.editCommandSet,'string',num2str(Comm.Set))
set(handles.editMeasureOnOff,'string',num2str(Measure.OnOff))
set(handles.editMeasureInput,'string',num2str(Measure.Input))
set(handles.editMeasurePrimary,'string',num2str(Measure.Electric))
set(handles.editMeasureSecondary,'string',num2str(Measure.Thermal))


function editCommandOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCommandSet_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandSet_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureSecondary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureSecondary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasurePrimary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasurePrimary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureInput_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureInput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)

function editMeasureOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function recordComm(handles)
global Plant
nG = length(Plant.Generator);
click = get(handles.CurrentComm,'String');

Comm.OnOff = str2double(get(handles.editCommandOnOff,'String'));
Com.Set = str2double(get(handles.editCommandSet,'String'));

Measure.OnOff = str2double(get(handles.editMeasureOnOff,'String'));
Measure.Input = str2double(get(handles.editMeasureInput,'String'));
Measure.Primary = str2double(get(handles.editMeasurePrimary,'String'));
Measure.Secondary = str2double(get(handles.editMeasureSecondary,'String'));
for j = 1:nG
    if strcmp(click, Plant.Generator(j).Name)%May need to change if multiple generators have the same name
        Plant.Generator(j).VariableStruct.Comm = Comm;
        Plant.Generator(j).VariableStruct.Measure = Measure;
    end
end

function popupmenuColor_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(hObject,'String');
color = list{get(hObject,'Value')};
colormap(handles.ResultPlot1,color);
colorVec = colormap(handles.ResultPlot1);
set(handles.figure1,'UserData',colorVec);

list = get(handles.uipanelMain1,'UserData');
nG = length(list);
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG));
for i=1:1:nG 
    set(handles.(strcat('Generator1_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Generator3_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
    set(handles.(strcat('Generator5_',num2str(i))),'BackgroundColor',colorsPlot(i,:))
end

function popupmenuColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
