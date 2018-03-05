function varargout = MainScreen1(varargin)
% MAINSCREEN1 MATLAB code for MainScreen1.fig
%      MAINSCREEN1, by itself, creates a new MAINSCREEN1 or raises the existing
%      singleton*.
%
%      H = MAINSCREEN1 returns the handle to a new MAINSCREEN1 or the handle to
%      the existing singleton*.
%
%      MAINSCREEN1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINSCREEN1.M with the given input arguments.
%
%      MAINSCREEN1('Property','Value',...) creates a new MAINSCREEN1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MainScreen1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MainScreen1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MainScreen1

% Last Modified by GUIDE v2.5 09-Feb-2018 10:01:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainScreen1_OpeningFcn, ...
                   'gui_OutputFcn',  @MainScreen1_OutputFcn, ...
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

% --- Executes just before MainScreen1 is made visible.
function MainScreen1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainScreen1 (see VARARGIN)

% Choose default command line output for MainScreen1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainScreen1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global Plant Model_dir SYSINDEX GENINDEX testSystems mainFig
mainFig = gcf;
movegui(gcf,'center');
set(gcf,'Name','Energy Planning Tool 2018.0.1')
if ~isfield(Plant,'Costs') 
    Costs = [];
else
    Costs = Plant.Costs;
end
Plant.Costs = defaultCosts(Costs,Plant.Generator);
if isempty(testSystems)
    SYSINDEX = 1;%index in the list of systems (full Plant structures)
else
    Names = {};
    for i = 1:1:length(testSystems)
        Names(end+1) = {testSystems(i).Name};
    end
    SYSINDEX = find(strcmp(Plant.Name,Names));
    if isempty(SYSINDEX)
        SYSINDEX = length(testSystems)+1;
    end
end
F = fieldnames(Plant);
for j = 1:1:length(F)
    testSystems(SYSINDEX).(F{j}) = Plant.(F{j});
end
GENINDEX = 1; %index in the list of generators of the current plant

%% Main Tabs
%Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2', etc.
TabText = {'Main Window';'Building Spec';'System Spec';'Network Spec';'Equip Costs';'Settings';};
for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        set(handles.uipanelMain1,'Visible','on')
        pos = get(handles.MainTab1','Position');
        set(handles.MainTab1,'Position',[pos(1),pos(2),pos(3),pos(4)+.5])
    else
        set(handles.(strcat('uipanelMain',j)),'Position',pan1pos,'Visible','off')
    end
end
set(handles.ProjectsPanel,'Position',[110, 42,100,6]);
handles = PlotAxes_Setup(hObject, eventdata, handles);
colormap(handles.axesMain,'parula');
set(hObject,'UserData',colormap(handles.axesMain));
set(handles.DesignDay,'value', 1);
%%setup library: give the buttons symbols
buttonIms = {'mGT', 'ICE', 'SOFC', 'MCFC', 'PEM', 'SolarPV', 'SolarThermal', 'Chiller', ...
            'AbChiller','AirHeater','WaterHeater', 'ColdStor', 'HotStor', 'Battery'}; % need images for: 'SolarStirling', 'Wind','HighTempStor
for i = 1:1:length(buttonIms)
    x = imread(fullfile(Model_dir,'GUI','Graphics',strcat(buttonIms{i}, '.png')));
    if license('test','Image Processing Toolbox')
        x = imresize(x,[42 70]);
    end
    set(handles.(buttonIms{i}), 'cdata', x)
end
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list)
handles = SysList_Make(handles);
popupmenu_Callback
set(handles.popupmenuBaseline,'String',list,'Value',1);


function handles = PlotAxes_Setup(hObject, eventdata, handles)
%primary axes
handles.axesMain = axes('Units','characters',...
    'Position', [13,16,125,22],'NextPlot','add',...
    'Tag', 'axesMain',...
    'Parent', handles.uipanelMain1,...
    'Visible','on');
%histogram axes
handles.axesCumulative = axes('Units','characters',...
    'Position', [155,16,50,22],'NextPlot','add','color','none',...
    'Tag', 'axesCumulative',...
    'Parent', handles.uipanelMain1,...
    'Visible','on');
%primary axes
handles.axesBuildLoad = axes('Units','characters',...
    'Position', [55,16,100,22],'NextPlot','add',...
    'Tag', 'axesBuildLoad',...
    'Parent', handles.uipanelMain2,...
    'Visible','on');
%histogram axes
handles.axesBuildCumulative = axes('Units','characters',...
    'Position', [170,16,40,22],'NextPlot','add','color','none',...
    'Tag', 'axesBuildCumulative',...
    'Parent', handles.uipanelMain2,...
    'Visible','on');


function popupmenu_Callback
tab = [];
handles = guihandles;
i = 1;
% Find out which tab is currently selected
while isempty(tab) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        tab = i;
    else
        i = i+1;
    end
end
set(handles.ProjectsPanel,'Visible','on')
if tab == 1
    loadMainTab(handles)
elseif tab == 2
    loadBuildingTab(handles)
elseif tab == 3
    Library_Callback([],[], handles)
    updateSystemRep(handles)
elseif tab == 4
    insertMockups(handles)
elseif tab == 5
    displayCosts
elseif tab == 6
    set(handles.ProjectsPanel,'Visible','off')
    loadSettingsTab(handles)
end

function loadMainTab(handles)
global TestData
if isfield(TestData,'Demand')
    list = {};
    if isfield(TestData.Demand,'E')
        list(end+1) = {'Electrical Demand'};
    end
    if isfield(TestData.Demand,'H')
        list(end+1) = {'Heating Demand'};
    end
    if isfield(TestData.Demand,'C')
        list(end+1) = {'Cooling Demand'};
    end
    if isfield(TestData.Demand,'W')
        list(end+1) = {'Water Demand'};
    end
    list(end+1) = {'Monthly Costs'};
elseif isfield(TestData,'Building')
    list = {'InternalGains';'NonHVACelectric';'Monthly Costs'};
end
set(handles.popupmenuAxes,'String',list,'Value',1);
days = max(2,floor(TestData.Timestamp(end) - TestData.Timestamp(1)));
if days<7
    set(handles.sliderZoom1,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data    
elseif days<31
    set(handles.sliderZoom1,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
elseif days<367
    set(handles.sliderZoom1,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
else
    set(handles.sliderZoom1,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
end
if days>20
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
else
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
end
set(handles.popupmenuAxes,'Value',1);
RunPlanning

function loadBuildingTab(handles)
global Model_dir testSystems SYSINDEX TestData
files = dir(fullfile(Model_dir, 'System Library','Buildings','*.mat'));
listB=strrep({files.name},'.mat','');
set(handles.popupmenuBuilding,'string',listB,'value',1);

files = dir(fullfile(Model_dir, 'System Library','Weather','*.mat'));
listW=strrep({files.name},'.mat','');
set(handles.popupmenuWeather,'string',listW,'value',1);
if isfield(testSystems(SYSINDEX),'Building') && ~isempty(testSystems(SYSINDEX).Building)
    I = find(strcmp(testSystems(SYSINDEX).Building.Name,listB));
    if isempty(I)
        listB(end+1) = {testSystems(SYSINDEX).Building.Name};
        I = length(listB);
    end
    set(handles.popupmenuBuilding,'string',listB);
    set(handles.popupmenuBuilding,'value',I);
    
    I = find(strcmp(TestData.Weather.Name,listW));
    if isempty(I)
        listW(end+1) = {TestData.Weather.Name};
        I = length(listW);
    end
    set(handles.popupmenuWeather,'string',listW);
    set(handles.popupmenuWeather,'value',I);
    
    set(handles.toggleHistorical,'Value',0);
    toggleSimulatedBuilding_Callback([], [], handles)
else
    toggleHistorical_Callback(handles);
end
days = max(2,floor(TestData.Timestamp(end) - TestData.Timestamp(1)));
if days<7   
    set(handles.sliderZoom2,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data   
elseif days<31
    set(handles.sliderZoom2,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
elseif days<367
    set(handles.sliderZoom2,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
else
    set(handles.sliderZoom2,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
end
if days>20
    set(handles.sliderDate2,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
else
    set(handles.sliderDate2,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
end
sliderDate2_Callback([],[],[])

function insertMockups(handles)
global Model_dir
handles.NetworkSetup = axes('Units','normalized',...
        'Position', [0,.15,1,.7],...
        'Tag', 'NetworkSetup',...
        'Parent', handles.uipanelMain4,...
        'Visible','off');
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','NetworkSetupMockup.png'));
image(x,'Parent',handles.NetworkSetup);
set(handles.NetworkSetup,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'box','off')

function displayCosts
global testSystems SYSINDEX GENINDEX
handles = guihandles;
quote='''';
nG = length(testSystems(SYSINDEX).Generator);
if ~isfield(testSystems(SYSINDEX),'Building')
    nB = 0;
else
    nB = length(testSystems(SYSINDEX).Building);
end
list = get(handles.uipanelMain5,'UserData');
if ~isempty(list)
    for i = 1:1:length(list(:,1))
        if isfield(handles,strcat('Equipment_',list{i,2}))
            delete(handles.(strcat('Equipment_',list{i,2})))
        end
    end
end
list = {};
for i = 1:1:(nG+nB)
    if i<=nG
        if ~strcmp(testSystems(SYSINDEX).Generator(i).Type,'Utility') && ~strcmp(testSystems(SYSINDEX).Generator(i).Type,'AC_DC')
            list(end+1,1:2) = {testSystems(SYSINDEX).Generator(i).Name, num2str(i)};
        end
    else
%         list(end+1,1:2) = {testSystems(SYSINDEX).Building(i-nG).Name, num2str(i)};
    end
end
set(handles.uipanelMain5,'UserData',list);
if length(list(:,1))>12
    set(handles.NextGen,'Visible','on','UserData',1)
    set(handles.PrevGen,'UserData',1,'Visible','off')
else
    set(handles.NextGen,'Visible','off')
    set(handles.PrevGen,'Visible','off')
end
%Makes buttons for GUI && status buttons next to corresponding generator
colorVec = get(gcf,'UserData');
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG+nB));
for i=1:1:length(list(:,1))
    curtab = floor((i-1)/12)+1;
    prev = 12*(curtab-1);
    if curtab==1
        vis = 'on';
    else
        vis = 'off';
    end
    pos = 45 - 2*(i-prev);
    k = str2double(list{i,2});
    callback = strcat('@(hObject,eventdata)MainScreen1(',quote,'SetGenCost_Callback',quote,',hObject,eventdata,guidata(hObject))');
    handles.(strcat('Equipment_',list{i,2})) = uicontrol('Style', 'pushbutton', 'String', list{i,1},...
    'Units','characters','Position', [1 pos 25 1.8],...
    'Tag', strcat('Equipment_',list{i,2}),'FontSize', 10,...
    'Parent', handles.uipanelMain5',...
    'Callback',eval(callback),'Visible',vis,'UserData',k,'BackgroundColor',colorsPlot(k,:));
end
set(handles.NPC_discount,'String',num2str(testSystems(SYSINDEX).Costs.DiscountRate));
updateSummaryTable(handles)
if GENINDEX>0
    SetGenCost_Callback([],[],[])
end

function loadSettingsTab(handles)
global SYSINDEX testSystems
set(handles.excessHeat, 'value', testSystems(SYSINDEX).optimoptions.excessHeat);
set(handles.excessCool, 'value', testSystems(SYSINDEX).optimoptions.excessCool);
set(handles.NoMixedInteger, 'value', ~testSystems(SYSINDEX).optimoptions.MixedInteger);
set(handles.MixedInteger, 'value', testSystems(SYSINDEX).optimoptions.MixedInteger);

set(handles.SpinReserve, 'value', testSystems(SYSINDEX).optimoptions.SpinReserve);
if testSystems(SYSINDEX).optimoptions.SpinReserve
    set(handles.SpinReservePerc, 'Visible','on','string', testSystems(SYSINDEX).optimoptions.SpinReservePerc);
else set(handles.SpinReservePerc, 'Visible','off');
end

nG = length(testSystems(SYSINDEX).Generator);
str = {};
stor = [];
for i = 1:1:nG
    if ismember(testSystems(SYSINDEX).Generator(i).Type,{'Electric Storage';'Thermal Storage'; 'Hydro Storage';})
        str{end+1} = testSystems(SYSINDEX).Generator(i).Name;
        stor(end+1) = i;
        if ~isfield(testSystems(SYSINDEX).Generator(i).VariableStruct,'Buffer')
            testSystems(SYSINDEX).Generator(i).VariableStruct.Buffer = 20;
        end
    end
end
if ~isempty(str)
    set(handles.StorageBuff, 'Visible', 'on','UserData',stor);
    set(handles.editBuffer, 'Visible', 'on');
    set(handles.textBuffer, 'Visible', 'on');
    set(handles.StorageBuff,'string',str,'value',1);
    set(handles.editBuffer, 'string', testSystems(SYSINDEX).Generator(stor(1)).VariableStruct.Buffer);
else
    set(handles.StorageBuff, 'Visible', 'off');
    set(handles.editBuffer, 'Visible', 'off');
    set(handles.textBuffer, 'Visible', 'off');
end

% --- Outputs from this function are returned to the command line.
function varargout = MainScreen1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
n = get(hObject,'Tag');
n = n(end);
m = [];
i = 1;
% Find out which tab is currently selected
while isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        m = i;
    else
        i = i+1;
    end
end
m = num2str(m);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% SaveBuilding color
bColor = get(handles.(strcat('MainTab',m)),'BackgroundColor');
set(handles.(strcat('MainTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% SaveBuilding dimensions
pos = get(handles.(strcat('MainTab',m)),'Position');
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.5])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.5])

%Change tab font
set(handles.(strcat('MainTab',m)),'FontWeight','normal','ForegroundColor',[0.501960784313726,0.501960784313726,0.501960784313726])
set(handles.(strcat('MainTab',n)),'FontWeight','bold','ForegroundColor',[0,0,0])

% SaveBuilding visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')
popupmenu_Callback

%% ProjectPanel Functions
function Save_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant Model_dir
Plant = testSystems(SYSINDEX);
[f,p]=uiputfile(fullfile(Model_dir,'Projects','New Project.mat'),'Save Project As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')

function Load_Callback(hObject, eventdata, handles)
global Plant Model_dir testSystems
cd(fullfile(Model_dir,'Projects'))
[fn,pn,~] = uigetfile('*.mat','Load Project File');
load(fullfile(pn,fn));
cd(Model_dir)
if ~isfield(Plant,'Costs') 
    Costs = [];
else
    Costs = Plant.Costs;
end
if isfield(testSystems,'Building') && ~isempty(testSystems(1).Building)
    Plant.Building = testSystems(1).Building;
end
Plant.Costs = defaultCosts(Costs,Plant.Generator);
newSys = length(testSystems)+1;
F = fieldnames(Plant);
for j = 1:1:length(F)
    testSystems(newSys).(F{j}) = Plant.(F{j});
end
handles = guihandles;
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list)
handles = SysList_Make(handles);
System_Callback(handles.(strcat('System_',num2str(newSys))), eventdata, handles)
popupmenu_Callback
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'));

function Copy_Callback(hObject, eventdata, handles)
global Plant testSystems SYSINDEX
Plant = testSystems(SYSINDEX);
Plant.Name = char(inputdlg('New Project Name','Select name',1,{strcat(Plant.Name,'_Alt')}));
testSystems(length(testSystems)+1) = Plant;
handles = guihandles;
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list);
handles = SysList_Make(handles);
System_Callback(handles.(strcat('System_',num2str(length(testSystems)))), eventdata, handles)
popupmenu_Callback
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'));

function System_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
if isfield(handles,strcat('System_',num2str(SYSINDEX)))
    pos = get(handles.(strcat('System_',num2str(SYSINDEX))),'Position');
    pos(4) = 1.8;
    set(handles.(strcat('System_',num2str(SYSINDEX))),'BackgroundColor',[1 1 .7],'ForegroundColor',[0.5 0.5 0.5],'Position',pos)
end
SYSINDEX = get(hObject,'UserData');
pos = get(handles.(strcat('System_',num2str(SYSINDEX))),'Position');
pos(4) = 2;
set(handles.(strcat('System_',num2str(SYSINDEX))),'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Position',pos)
pos2 = get(handles.Save,'Position');
pos2(1) = pos(1);
set(handles.Save,'Position',pos2);
GENINDEX = min(length(testSystems(SYSINDEX).Generator), GENINDEX); % prevent index exceeds matrix dimensions error
popupmenu_Callback

function DeleteSys_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
GENINDEX = 1;
SYSINDEX = get(hObject,'UserData');
val = get(handles.popupmenuBaseline,'Value');
if SYSINDEX<=val
    val = max(1,val-1);
end
if length(testSystems)>SYSINDEX
    testSystems(SYSINDEX:length(testSystems)-1) = testSystems(SYSINDEX+1:length(testSystems));
end
testSystems = testSystems(1:length(testSystems)-1);
handles = guihandles;
handles = SysList_Make(handles);
System_Callback(handles.(strcat('System_',num2str(1))), eventdata, handles)
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
set(handles.ProjectsPanel,'UserData',list);
set(handles.popupmenuAxes,'Value',1);
RunPlanning
set(handles.popupmenuBaseline,'String',get(handles.ProjectsPanel,'UserData'),'Value',val);

function PrevSys_Callback(hObject, eventdata, handles)
list = get(handles.ProjectsPanel,'UserData');
gen = length(list);
page = get(handles.PrevSys,'UserData');%current page of the list
if page<2
    set(handles.PrevSys,'Visible','off','UserData',page-1)
else
    set(handles.PrevSys,'Visible','on','UserData',page-1);
end
set(handles.NextSys,'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('System_',j)),'Visible','off');
         set(handles.(strcat('DeleteSys',j)),'Visible','off');
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('System_',j)),'Visible','on');
     set(handles.(strcat('DeleteSys',j)),'Visible','on');
end

function NextSys_Callback(hObject, eventdata, handles)
list = get(handles.ProjectsPanel,'UserData');
gen = length(list);
page = get(handles.PrevSys,'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.NextSys,'Visible','off','UserData',page+1)
else
    set(handles.NextSys,'Visible','on','UserData',page+1);
end
set(handles.PrevSys,'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('System_',j)),'Visible','off');
     set(handles.(strcat('DeleteSys',j)),'Visible','off');
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('System_',j)),'Visible','on');
         set(handles.(strcat('DeleteSys',j)),'Visible','on');
    end
end

%% Tab 1 functions
function pushbuttonEDC_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant mainFig
mainFig = [];  
Plant = testSystems(SYSINDEX);
close
DISPATCH

function popupmenuAxes_Callback(hObject, eventdata, handles)
RunPlanning

function PlotDemand(handles)
list = get(handles.popupmenuDemand,'String');
item = list{get(handles.popupmenuDemand,'Value')};
tab = [];
i = 1;
% Find out which tab is currently selected
while isempty(tab) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        tab = i;
    else
        i = i+1;
    end
end
PlotHistorical(handles,item,tab)

% --- Executes during object creation, after setting all properties.
function popupmenuAxes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderDate1_Callback(hObject, eventdata, handles)
RunPlanning

% --- Executes during object creation, after setting all properties.
function sliderDate1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function sliderZoom1_Callback(hObject, eventdata, handles)
if get(handles.sliderZoom1,'Value')== get(handles.sliderZoom1,'Max')
    set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
elseif strcmp(get(handles.sliderDate1,'Visible'),'off')
    set(handles.sliderDate1,'Visible','on');set(handles.textDate1,'Visible','on');
end
RunPlanning

% --- Executes during object creation, after setting all properties.
function sliderZoom1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Tab 2 functions (any changes to project on this tab require re-running design days)
function popupmenuBuilding_Callback(hObject, eventdata, handles)
%load the selected building
global testSystems SYSINDEX Model_dir
list = get(handles.popupmenuBuilding,'string');
sel = get(handles.popupmenuBuilding,'Value');
load(fullfile(Model_dir,'System Library','Buildings',list{sel}));
testSystems(SYSINDEX).Building = building;
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])
testSystems(SYSINDEX).Design = [];%empty design day solution

function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuWeather_Callback(hObject, eventdata, handles)
%load selected weather profile
global Model_dir testSystems TestData
list = get(handles.popupmenuWeather,'string');
sel = get(handles.popupmenuWeather,'Value');
load(fullfile(Model_dir,'System Library','Weather',list{sel}));
TestData.Weather = interpolateWeather(weather,TestData.Timestamp);
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])
for i_ts = 1:1:length(testSystems)
    testSystems(i_ts).Design = [];%empty design day solution for all cases
end

function popupmenuWeather_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in toggleSimulatedBuilding.
function toggleSimulatedBuilding_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Model_dir TestData
if ~isfield(TestData,'Weather') || isempty(TestData.Weather)
    list = get(handles.popupmenuWeather,'string');
    sel = get(handles.popupmenuWeather,'Value');
    load(fullfile(Model_dir,'System Library','Weather',list{sel}));
    loadTestWeather(weather);
    popupmenuBuilding_Callback([], [], handles)
end
list = {};
if isfield(TestData,'Demand')
    if isfield(TestData.Demand,'E')
        list(end+1,1) = {'Electrical Demand'};
    end
    if isfield(TestData.Demand,'H')
        list(end+1,1) = {'Heating Demand'};
    end
    if isfield(TestData.Demand,'C')
        list(end+1,1) = {'Cooling Demand'};
    end
    if isfield(TestData.Demand,'W')
        list(end+1,1) = {'Water Demand'};
    end
end
if isfield(TestData,'Building')
    list(end+1:end+2) = {'InternalGains';'NonHVACelectric';};
end
set(handles.popupmenuDemand,'String',list);
if get(handles.toggleHistorical,'Value')==1
    a = get(handles.toggleHistorical,'BackgroundColor');
    b = get(handles.toggleSimulatedBuilding,'BackgroundColor');
    c = get(handles.toggleHistorical,'ForegroundColor');
    d = get(handles.toggleSimulatedBuilding,'ForegroundColor');
    set(handles.toggleSimulatedBuilding,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.toggleHistorical,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.toggleSimulatedBuilding,'Value',1); %was already pressed
end
list = {'textBuilding';'pushbuttonPrev';'pushbuttonNext';'popupmenuBuilding';'popupmenuWeather';'SaveBuilding';'pushbuttonReset';...
        'textVariability';'textDaily';'editDailyVariability';'textInterHour';'editHourlyVariability';'textIntraHour';'editSubHourlyVariability';};
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','On');
end
set(handles.LoadTestData,'Visible','Off')
pushbuttonPrev_Callback([], [], handles);
set(handles.editDailyVariability,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.DailyVariability))
set(handles.editHourlyVariability,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.HourlyVariability))
set(handles.editSubHourlyVariability,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.SubHourlyVariability))

% --- Executes on button press in toggleHistorical.
function toggleHistorical_Callback(handles)
global TestData
Outs =  fieldnames(TestData.Demand);
list = {};
if strcmp(Outs{1},'E')
    list(end+1) = {'Electrical Demand'};
elseif strcmp(Outs{1},'C')
    list(end+1) = {'Cooling Demand'};
elseif strcmp(Outs{1},'H')
    list(end+1) = {'Heating Demand'};
elseif strcmp(Outs{1},'W')
    list(end+1) = {'Water Demand'};
end
set(handles.popupmenuDemand,'String',list,'Value',1);
list(end+1) = {'Monthly Costs'};
set(handles.popupmenuAxes,'String',list,'Value',1);
if get(handles.toggleSimulatedBuilding,'Value')==1
    a = get(handles.toggleHistorical,'BackgroundColor');
    b = get(handles.toggleSimulatedBuilding,'BackgroundColor');
    c = get(handles.toggleHistorical,'ForegroundColor');
    d = get(handles.toggleSimulatedBuilding,'ForegroundColor');
    set(handles.toggleSimulatedBuilding,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.toggleHistorical,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else
    set(handles.toggleHistorical,'Value',1); %was already pressed
end
set(handles.uipanelBuildParam2,'Visible','Off');
set(handles.uipanelBuildParam1,'Visible','Off');
list = {'textBuilding';'pushbuttonPrev';'pushbuttonNext';'popupmenuBuilding';'popupmenuWeather';'SaveBuilding';'pushbuttonReset';...
        'textVariability';'textDaily';'editDailyVariability';'textInterHour';'editHourlyVariability';'textIntraHour';'editSubHourlyVariability';};
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','Off');
end
set(handles.LoadTestData,'Visible','On')


function pushbuttonPrev_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
set(handles.uipanelBuildParam2,'Visible','Off');
set(handles.uipanelBuildParam1,'Visible','On');
set(handles.editArea,'String',num2str(testSystems(SYSINDEX).Building.Area*10.76))
set(handles.editWindowWall,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.WindowWallRatio))
set(handles.editOccupancy,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.occupancy*10.76))
set(handles.editLighting,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.InteriorLights*10.76))
set(handles.editEquipment,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.equipment*10.76))

function pushbuttonNext_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
set(handles.uipanelBuildParam2,'Visible','On');
set(handles.uipanelBuildParam1,'Visible','Off');
set(handles.editComfort,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.Comfort*9/5))
set(handles.editRvalue,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.WallRvalue*1055/3600*9/5*10.76/1000))
set(handles.editUvalue,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.WindowUvalue))
set(handles.editAirChange,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.AirChangePerHr))
set(handles.editDewPoint,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.DPset*9/5+32))
set(handles.editColdSupply,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.ColdAirSet*9/5+32))
set(handles.editFanPower,'String',num2str(testSystems(SYSINDEX).Building.VariableStruct.FanPower))


% --- Executes on button press in LoadTestData.
function LoadTestData_Callback
global TestData testSystems
TestData = [];
%load historical data
loadTestData
for i_ts = 1:1:length(testSystems)
    testSystems(i_ts).Design = [];%empty design day solution
end


function SaveBuilding_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Model_dir
%savebuilding building type
[f,p]=uiputfile(fullfile(Model_dir,'System Library','Buildings',strcat(testSystems(SYSINDEX).Building.Name,'.mat')),'Save Building As...');
if f==0; return; end
testSystems(SYSINDEX).Building.Name = strrep(f,'.mat','');
building = testSystems(SYSINDEX).Building;
save([p,f],'building')

function pushbuttonReset_Callback(hObject, eventdata, handles)
popupmenuBuilding_Callback(hObject, eventdata, handles)
pushbuttonPrev_Callback(hObject, eventdata, handles);


function editArea_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.Area = str2double(get(hObject,'String'))/10.76;%convert to m^2
% volume of treated air space (assume height of 3m)
testSystems(SYSINDEX).Building.Volume = testSystems(SYSINDEX).Building.Area*3; % m^3
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editArea_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonOccupancy_Callback(hObject, eventdata, handles)
global testSystems
handles = guihandles;
ScheduleEditor(hObject,'occ',handles)
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],handles)

function editOccupancy_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.occupancy = str2double(get(hObject,'String'))/10.76;%convert to m^2
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editOccupancy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonLighting_Callback(hObject, eventdata, handles)
global testSystems
handles = guihandles;
ScheduleEditor(hObject,'light',handles)
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],handles)

function editLighting_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.lighting = str2double(get(hObject,'String'))/10.76;%convert to m^2
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editLighting_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonEquipment_Callback(hObject, eventdata, handles)
global testSystems
handles = guihandles;
ScheduleEditor(hObject,'equip',handles)
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],handles)

function editEquipment_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.equipment = str2double(get(hObject,'String'))/10.76;%convert to m^2
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editEquipment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editComfort_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.Comfort = str2double(get(hObject,'String'))*5/9;%convert to Celcius
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editComfort_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editRvalue_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
%aproximate R-values: Windows = 3, doors = 5, walls = .03, floor/ceiling 50
% U = 1/R-value
% Q= U*A*deltaT, resistance = 1/U
% a value of U = 0.2 BTU/hr-F-ft^2 converts to 1.1352e-3 kJ/s-K-m^2 (0.2*1055/3600*9/5*10.76/1000)  10.76ft^2 per m^2, 1055J/BTU, 9/5F per K
% inverting gives Resistance = 880.92 m^2*K/kW
testSystems(SYSINDEX).Building.VariableStruct.WallRvalue = str2double(get(hObject,'String'))*3600/1055*5/9*1000/10.76;% convert hr-F-ft^2/BTU to m^2*K/kW
testSystems(SYSINDEX).Building.VariableStruct.RoofRvalue = str2double(get(hObject,'String'))*3600/1055*5/9*1000/10.76;% convert hr-F-ft^2/BTU to m^2*K/kW
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editRvalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAirChange_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.AirChangePerHr = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editAirChange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDewPoint_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.DPset = (str2double(get(hObject,'String'))-32)*5/9;
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editDewPoint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editColdSupply_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.ColdAirSet = (str2double(get(hObject,'String'))-32)*5/9;
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editColdSupply_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editFanPower_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.Fan_Power = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editFanPower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWindowWall_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.WindowWallRatio = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editWindowWall_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editUvalue_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.WindowUvalue = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editUvalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDailyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.DailyVariability = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editDailyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editHourlyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.HourlyVariability = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editHourlyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSubHourlyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.SubHourlyVariability = str2double(get(hObject,'String'));
reLoadBuilding(testSystems(1).Building,testSystems(1).Network)
sliderDate2_Callback([],[],[])

function editSubHourlyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection savebuilding in popupmenuDemand.
function popupmenuDemand_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.popupmenuDemand,'String');
item = list{get(handles.popupmenuDemand,'Value')};
PlotHistorical(handles,item,2)

% --- Executes during object creation, after setting all properties.
function popupmenuDemand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderDate2_Callback(hObject, eventdata, handles)
if isempty(handles) || strcmp(get(hObject,'Tag'),'sliderDate2')
    handles = guihandles;
end
PlotDemand(handles)

function sliderDate2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderZoom2_Callback(hObject, eventdata, handles)
handles = guihandles;
if get(handles.sliderZoom2,'Value')== get(handles.sliderZoom2,'Max')
    set(handles.sliderDate2,'Visible','off');set(handles.textDate2,'Visible','off');
elseif strcmp(get(handles.sliderDate2,'Visible'),'off')
    set(handles.sliderDate2,'Visible','on');set(handles.textDate2,'Visible','on');
end
sliderDate2_Callback([],[],[])

function sliderZoom2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% Tab 3 functions (any changes to the system specification or perfomance features requires the designn days to be re-run)
function saveSystem_Callback(hObject, eventdata, handles)
global GENINDEX SYSINDEX testSystems Model_dir
component = testSystems(SYSINDEX).Generator(GENINDEX);
[f,p]=uiputfile(fullfile(Model_dir,'System Library',component.Name),'Save Plant As...');
if f==0; return; end
save([p,f],'component')

% --- Executes on button press in Library.
function Library_Callback(hObject, eventdata, handles)
set(handles.uipanelLibrary,'Visible','on');
set(handles.uipanelGenSpec,'Visible','off');
set(handles.Library,'Visible','off');
set(handles.saveSystem,'Visible','off');
set(handles.pushbuttonRemove,'Visible','off');

%% The following need to represent different things for each system selected in 3rd tab
function CompName_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).Name = get(hObject,'String');
old = char(testSystems(SYSINDEX).Network.Equipment(GENINDEX));
first = strfind(old,'.');
old = old(1:first);
new = cellstr(strcat(old,get(hObject,'String')));
testSystems(SYSINDEX).Network.Equipment(GENINDEX) = new;

function CompName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function compText1_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(1,1) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'AC_DC')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.AC_to_DC_eff = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'}) 
    handles = guihandles;
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'UB',str2double(get(hObject,'String')));
    p = eig(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.StateSpace.A);
    w_0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w_0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    secondOrderResponse(testSystems(SYSINDEX).Generator(GENINDEX),handles);
elseif ismember(Gen.Type,{'Electric Storage';'Thermal Storage';}) 
    testSystems(SYSINDEX).Generator(GENINDEX).Size = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Size = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText2_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(1,2) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'AC_DC')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.DC_to_AC_eff = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'}) 
    handles = guihandles;
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'LB',str2double(get(hObject,'String')));
    p = eig(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    secondOrderResponse(testSystems(SYSINDEX).Generator(GENINDEX),handles);
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Sizem2 = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Voltage = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SizeLiter = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText3_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(2,1) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'AC_DC')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Capacity = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller'})
    handles = guihandles;
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'dX_dt',str2double(get(hObject,'String')));
    p = eig(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    set(handles.compText6,'String',num2str(zeta));
    secondOrderResponse(testSystems(SYSINDEX).Generator(GENINDEX),handles);
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Eff = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MaxDOD = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tcold = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText4_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(2,2) = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.StartCost = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Azimuth = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.ChargeResist = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Thot = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText5_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(3,1) = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    handles = guihandles;
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'w0',str2double(get(hObject,'String')));
    set(handles.compText3,'String',num2str(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.dX_dt));
    secondOrderResponse(testSystems(SYSINDEX).Generator(GENINDEX),handles);
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tilt = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.DischResist = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.dX_dt = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText6_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(3,2) = str2double(get(hObject,'String'));
elseif ismember(Gen.Type,{'CHP Generator';'Electric Generator';'Heater';'Chiller';})
    handles = guihandles;
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'zeta',str2double(get(hObject,'String')));
    p = eig(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.StateSpace.A);
    w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
    zeta = -real(p(1)+p(2))/(2*w0);
    set(handles.compText5,'String',num2str(w0));
    secondOrderResponse(testSystems(SYSINDEX).Generator(GENINDEX),handles);
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PeakCharge = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.ChargeEff = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText7_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PeakDisch = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.DischargeEff = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText8_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SelfDischarge = str2double(get(hObject,'String'))/(31*24*100);
elseif strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SelfDischarge = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText9_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.FillRate = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function compText10_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Thermal Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.DischRate = str2double(get(hObject,'String'));
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function uitableEffCurve_CellEditCallback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Outputs = fieldnames(testSystems(SYSINDEX).Generator(GENINDEX).Output);
handlesC=get(handles.uipanelGenSpec,'Children');
nOutput = eventdata.Indices;
newValue = eventdata.NewData;
if length(Outputs) == 5
    for i= 1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'checkboxSeasonal')
            if get(handlesC(i),'Value') == 1
                type = get(handlesC(i-1),'Value');
                if type == 1
                    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRateTable(nOutput(1),nOutput(2)) = newValue; 
                else
                    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.WinRateTable(nOutput(1),nOutput(2)) = newValue;
                end
            else
                testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRateTable(nOutput(1),nOutput(2)) = newValue;
            end
        end
    end
else
    testSystems(SYSINDEX).Generator(GENINDEX).Output.(Outputs{nOutput(2)})(nOutput(1)) = newValue;
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function popupRates_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
handlesC = get(handles.uipanelGenSpec,'Children');

if hObject.Value == 1
    for i=1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'compText1')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(1,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText2')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(1,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText3')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(2,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText4')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(2,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText5')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(3,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText6')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.SumRates(3,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'uitableEffCurve')
            set(handlesC(i),'Data',Gen.VariableStruct.SumRateTable);
        end
    end
    Start1 = Gen.VariableStruct.SumStartMonth;
    Start2 = Gen.VariableStruct.SumStartDay;
    if Gen.VariableStruct.WinStartDay == 1
        End1 = Gen.VariableStruct.WinStartMonth-1;
        End2 = m_d(End1,2);
    else
        End1 = Gen.VariableStruct.WinStartMonth;
        End2 = Gen.VariableStruct.WinStartDay - 1;
    end
    tag1 = 'popupDates1S';
    tag2 = 'popupDates2S';
    tag3 = 'popupDates3S';
    tag4 = 'popupDates4S';
else
    for i=1:length(handlesC)
        if strcmp(get(handlesC(i),'Tag'),'compText1')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(1,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText2')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(1,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText3')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(2,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText4')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(2,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText5')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(3,1)));
        elseif strcmp(get(handlesC(i),'Tag'),'compText6')
            set(handlesC(i),'String',num2str(Gen.VariableStruct.WinRates(3,2)));
        elseif strcmp(get(handlesC(i),'Tag'),'uitableEffCurve')
            set(handlesC(i),'Data',Gen.VariableStruct.WinRateTable);
        end
    end
    Start1 = Gen.VariableStruct.WinStartMonth;
    Start2 = Gen.VariableStruct.WinStartDay;
    if Gen.VariableStruct.SumStartDay == 1
        End1 = Gen.VariableStruct.SumStartMonth-1;
        End2 = m_d(End1,2);
    else
        End1 = Gen.VariableStruct.SumStartMonth;
        End2 = Gen.VariableStruct.SumStartDay - 1;
    end
    tag1 = 'popupDates1W';
    tag2 = 'popupDates2W';
    tag3 = 'popupDates3W';
    tag4 = 'popupDates4W';
end
days=(1:m_d(Start1,2))';
EditSystem('createPopup',handles,{'1','2','3','4','5','6','7','8','9','10','11','12'},...
            [10 29 7 1],12,'normal',tag1,Start1,'popupDates')
EditSystem('createPopup',handles,days,[18 29 7 1],12,'normal',tag2,Start2,'popupDates')
clear days
days=(1:m_d(End1,2))';
EditSystem('createPopup',handles,{'1','2','3','4','5','6','7','8','9','10','11','12'},...
            [30 29 7 1],12,'normal',tag3,End1,'popupDates')
EditSystem('createPopup',handles,days,[38 29 7 1],12,'normal',tag4,End2,'popupDates')
testSystems(SYSINDEX).Design = [];%empty design day solution

function popupDates_Callback(hObject,handles0,handles)
global GENINDEX SYSINDEX testSystems
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
handlesc = get(handles.uipanelGenSpec,'Children');
xW=[];
xS=[];

for i=1:(length(handlesc))
    if ~isempty(strfind(handlesc(i).Tag,'popupDates')) && ~isempty(strfind(handlesc(i).Tag,'S'))
        xS(:,end+1)=[get(handlesc(i),'Value');i];
    elseif ~isempty(strfind(handlesc(i).Tag,'popupDates')) && ~isempty(strfind(handlesc(i).Tag,'W'))
        xW(:,end+1)=[get(handlesc(i),'Value');i];
    end
end
Date = get(hObject,'Value');
Type = get(hObject,'Tag');
m_d=[1 31; 2 28; 3 31; 4 30; 5 31; 6 30; 7 31; 8 31; 9 30; 10 31; 11 30; 12 31];


if ismember(Type,{'popupDates1S','popupDates2S','popupDates3S','popupDates4S'})
    xS = xS(:,1:4);
    xS = fliplr(xS);
    if strcmp(Type,'popupDates1S')
        Gen.VariableStruct.SumStartMonth = Date;
    elseif strcmp(Type,'popupDates2S')
        Gen.VariableStruct.SumStartDay = Date;
    elseif strcmp(Type,'popupDates3S')
        if xS(1,4) >= m_d(xS(1,3),2)
            Gen.VariableStruct.WinStartMonth = Date + 1;
            Gen.VariableStruct.WinStartDay = 1;
            if Gen.VariableStruct.WinStartMonth == 13
                Gen.VariableStruct.WinStartMonth = 1;
            end
        else
            Gen.VariableStruct.WinStartMonth = Date;
            Gen.VariableStruct.WinStartDay = xS(1,4)+1;
        end
    elseif strcmp(Type,'popupDates4S')
        if m_d(xS(1,3),2) == Date
            Gen.VariableStruct.WinStartMonth = xS(1,3) + 1;
            Gen.VariableStruct.WinStartDay = 1;
            if Gen.VariableStruct.WinStartMonth == 13
                Gen.VariableStruct.WinStartMonth = 1;
            end
        else
            Gen.VariableStruct.WinStartMonth = xS(1,3);
            Gen.VariableStruct.WinStartDay = Date + 1;
        end
    end
    days=(1:m_d(xS(1,1),2))';
    if get(handlesc(xS(2,2)),'Value')>length(days)
        set(handlesc(xS(2,2)),'Value',days(end))
    end
    set(handlesc(xS(2,2)),'String',days);
    clear days
    days=(1:m_d(xS(1,3),2))';
    if get(handlesc(xS(2,4)),'Value')>length(days)
        set(handlesc(xS(2,4)),'Value',days(end))
    end
    set(handlesc(xS(2,4)),'String',days);
elseif ismember(Type,{'popupDates1W','popupDates2W','popupDates3W','popupDates4W'})
    xW = xW(:,1:4);
    xW = fliplr(xW);
    if strcmp(Type,'popupDates1W')
        Gen.VariableStruct.WinStartMonth = Date;
    elseif strcmp(Type,'popupDates2W')
        Gen.VariableStruct.WinStartDay = Date;
    elseif strcmp(Type,'popupDates3W')
        if xW(1,4) >= m_d(xW(1,3),2)
            Gen.VariableStruct.SumStartMonth = Date + 1;
            Gen.VariableStruct.SumStartDay = 1;
            if Gen.VariableStruct.SumStartMonth == 13
                Gen.VariableStruct.SumStartMonth = 1;
            end
        else
            Gen.VariableStruct.SumStartMonth = Date;
            Gen.VariableStruct.SumStartDay = xW(1,4)+1;
        end
    elseif strcmp(Type,'popupDates4W')
        if m_d(xW(1,3),2) == Date
            Gen.VariableStruct.SumStartMonth = xW(1,3) + 1;
            Gen.VariableStruct.SumStartDay = 1;
            if Gen.VariableStruct.SumStartMonth == 13
                Gen.VariableStruct.SumStartMonth = 1;
            end
        else
            Gen.VariableStruct.SumStartMonth = xW(1,3);
            Gen.VariableStruct.SumStartDay = Date + 1;
        end
    end
    days=(1:m_d(xW(1,1),2))';
    if get(handlesc(xW(2,2)),'Value')>length(days)
        set(handlesc(xW(2,2)),'Value',days(end))
    end
    set(handlesc(xW(2,2)),'String',days);
    clear days
    days=(1:m_d(xW(1,3),2))';
    if get(handlesc(xW(2,4)),'Value')>length(days)
        set(handlesc(xW(2,4)),'Value',days(end))
    end
    set(handlesc(xW(2,4)),'String',days);
end
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct=Gen.VariableStruct;
testSystems(SYSINDEX).Design = [];%empty design day solution
              
function radiobuttonGridsellback_Callback(hObject,handles0,handles)
global GENINDEX SYSINDEX testSystems
handles = guihandles;
if strcmp(get(hObject,'Tag'),'None')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = 0;
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc = 0;
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh = 0;
    set(handles.textEdit11,'String','Minimum Import (kW)');
    set(handles.SellbackAsPerc,'Value',0);
    set(handles.FixedRate,'Value',0);
    set(handles.maxSellback,'String','0');
elseif strcmp(get(hObject,'Tag'),'SellbackAsPerc')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = -1;
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc = str2double(get(handles.editTariffs,'String'));
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh = -inf;
    set(handles.textEdit11,'String','Maximum Export (kW)');
    set(handles.None,'Value',0);
    set(handles.FixedRate,'Value',0);
    set(handles.maxSellback,'String','inf');
elseif strcmp(get(hObject,'Tag'),'FixedRate')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = str2double(get(handles.editSellbackRate,'String'));
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc = 0;
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh = -inf;
    set(handles.textEdit11,'String','Maximum Export (kW)');
    set(handles.SellbackAsPerc,'Value',0);
    set(handles.None,'Value',0);
    set(handles.maxSellback,'String','inf');
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function editTariffs(hObject,handles)
global GENINDEX SYSINDEX testSystems
handles = guihandles;
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = -1;
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc = str2double(get(handles.editTariffs,'String'));

function editSellbackRate(hObject,handles)
global GENINDEX SYSINDEX testSystems
handles = guihandles;
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = str2double(get(handles.editSellbackRate,'String'));
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc = 0;

function maxSellback(hObject,handles)
global GENINDEX SYSINDEX testSystems
handles = guihandles;
if get(handles.None,'Value') == 0
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh = - str2double(get(handles.editMaxSellback,'String'));
else
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh = str2double(get(handles.editMaxSellback,'String'));
end


function radioSource_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
handlesc = get(handles.uipanelGenSpec,'Children');
for i =1:length(handlesc)
    if strcmp(get(handlesc(i),'Tag'),'radioSource')
        if strcmp(get(handlesc(i),'String'),get(hObject,'String'))
            set(handlesc(i),'Value',1)
        else
            set(handlesc(i),'Value',0)
        end
    end
end
if strcmp(get(hObject,'String'),'Natural Gas')
    testSystems(SYSINDEX).Generator(GENINDEX).Source = 'NG';
else
    testSystems(SYSINDEX).Generator(GENINDEX).Source = get(hObject,'String');
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function popupSolar_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = get(hObject,'Value');
state = stateName(stateNum);
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.State = char(state);
testSystems(SYSINDEX).Design = [];%empty design day solution

function radiobuttonSType_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
handlesc = get(handles.uipanelGenSpec,'Children');
for i =1:length(handlesc)
    if strcmp(get(handlesc(i),'Tag'),'radiobuttonSType')
        if strcmp(get(handlesc(i),'String'),get(hObject,'String'))
            set(handlesc(i),'Value',1)
        else
            set(handlesc(i),'Value',0)
        end
    end
end
if strcmp(get(hObject,'String'),'Flat Panel')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PVtype = 'flat';
elseif strcmp(get(hObject,'String'),'Concentrated')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PVtype = 'concentrated';
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function radiobuttonSTracking_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
handlesc = get(handles.uipanelGenSpec,'Children');
for i =1:length(handlesc)
    if strcmp(get(handlesc(i),'Tag'),'radiobuttonSTracking')
        if strcmp(get(handlesc(i),'String'),get(hObject,'String'))
            set(handlesc(i),'Value',1)
        else
            set(handlesc(i),'Value',0)
        end
    end
end
if strcmp(get(hObject,'String'),'Fixed')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tracking = 'fixed';
elseif strcmp(get(hObject,'String'),'Single Axis')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tracking = '1axis';
elseif strcmp(get(hObject,'String'),'Dual Axis')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tracking = '2axis';
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function uitableDCAC_CellEditCallback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Data = get(hObject,'Data');
testSystems(SYSINDEX).Design = [];%empty design day solution

function uitableBat_CellEditCallback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.VoltCurve = get(hObject,'Data');
testSystems(SYSINDEX).Design = [];%empty design day solution

function AddSystem(hObject, eventdata,handles)
%AddSystem(hObject,handles)
%Adds a generator or storage system to the plant.
global testSystems SYSINDEX Model_dir GENINDEX
nG = length(testSystems(SYSINDEX).Generator);
GENINDEX = nG+1;
emptyStr = '[None]';
files=dir(fullfile(Model_dir,'System Library',get(hObject,'Tag'),'*.mat'));
list=strrep({files.name},'.mat','');
if isempty(list)
    list{1} = emptyStr;
end
[s,OK] = listdlg('PromptString','Select Model', ...
    'SelectionMode','single', ...
    'ListString',list);
if OK && ~strcmp(list{s},emptyStr)
    componentName = list{s};
    load(fullfile(Model_dir,'System Library',get(hObject,'Tag'), ...
        strcat(componentName,'.mat')));
    if strcmp(component.Type,'Utility')
        testSystems(SYSINDEX).Generator(GENINDEX) = struct( ...
            'Type',component.Type, ...
            'Name','Elec Utility1', ...
            'Source','Electricity', ...
            'Output',struct('Capacity',[], ...
                'Electricity',[], ...
                'Heat',[], ...
                'Steam',[], ...
                'Cooling',[]), ...
            'Size',1, ...
            'Enabled',1, ...
            'VariableStruct',rmfield(component, {'Type','Source'}));
    else
        if isfield(testSystems(SYSINDEX).Generator(1),'QPform')
            component.QPform = [];
        end
        testSystems(SYSINDEX).Generator(GENINDEX) = component;
    end
    testSystems(SYSINDEX).Network(1).Equipment(end+1) = {strcat( ...
    testSystems(SYSINDEX).Generator(GENINDEX).Type,'.',testSystems(SYSINDEX).Generator(GENINDEX).Name)};
    EditSystem(handles)
end
updateSystemRep(handles)
updateSystemList
if ~isfield(testSystems(SYSINDEX),'Costs') 
    Costs = [];
else
    Costs = testSystems(SYSINDEX).Costs;
end
testSystems(SYSINDEX).Costs = defaultCosts(Costs,testSystems(SYSINDEX).Generator);
testSystems(SYSINDEX).Design = [];%empty design day solution

% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX 
if GENINDEX > 0
    str = strcat(testSystems(SYSINDEX).Generator(GENINDEX).Type,'.',testSystems(SYSINDEX).Generator(GENINDEX).Name);
    for n = 1:1:length(testSystems(SYSINDEX).Network)
        testSystems(SYSINDEX).Network(n).Equipment = testSystems(SYSINDEX).Network(n).Equipment(~strcmp(str,testSystems(SYSINDEX).Network(n).Equipment));
    end
    testSystems(SYSINDEX).Generator = testSystems(SYSINDEX).Generator([1:GENINDEX-1,GENINDEX+1:length(testSystems(SYSINDEX).Generator)]);
    testSystems(SYSINDEX).Costs.Equipment = testSystems(SYSINDEX).Costs.Equipment([1:GENINDEX-1,GENINDEX+1:length(testSystems(SYSINDEX).Costs.Equipment)]);
    GENINDEX = 1;
    updateSystemRep(handles)
    updateSystemList
    Library_Callback([],[], handles);
end
testSystems(SYSINDEX).Design = [];%empty design day solution

function CompFuel_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
fuels = get(hObject,'String');
testSystems(SYSINDEX).Generator(GENINDEX).Source = fuels{get(hObject,'Value')};
testSystems(SYSINDEX).Design = [];%empty design day solution

function CompFuel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tab 5 functions
function PrevGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain5,'UserData');
gen = length(list);
page = get(handles.PrevGen,'UserData');%current page of the list
if (page-1)<2%if the new page is the 1st
    set(handles.PrevGen,'Visible','off','UserData',page-1)
else
    set(handles.PrevGen,'Visible','on','UserData',page-1);
end
set(handles.NextGen,'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('Equipment_',j)),'Visible','off');
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('Equipment_',j)),'Visible','on');
end

function NextGen_Callback(hObject, eventdata, handles)
handles = guihandles;
list = get(handles.uipanelMain5,'UserData');
gen = length(list);
page = get(handles.PrevGen,'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.NextGen,'Visible','off','UserData',page+1)
else
    set(handles.NextGen,'Visible','on','UserData',page+1);
end
set(handles.PrevGen,'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('Equipment_',j)),'Visible','off');
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('Equipment_',j)),'Visible','on');
    end
end

function SetGenCost_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
if ~isempty(hObject)
    GENINDEX = get(hObject,'UserData');
end
handles = guihandles;
list = {'Name';'Size';'Cost';'CostkW';'O_M';'Financed';'Loan';'LoanTerm';};
for i = 1:1:length(list(:,1))
    if isfield(handles,strcat('Cost_',list{i},'_text'))
        delete(handles.(strcat('Cost_',list{i},'_text')))
        delete(handles.(strcat('Cost_',list{i},'_edit')))
        delete(handles.(strcat('Cost_',list{i},'_unit')))
    end
end
if isfield(handles,'OptimizeGen')
    delete(handles.OptimizeGen);
end
listText = {'Name';'Size';'Total Cost';'Normalized Cost';'Operations & Maintenance';'Percent Financed';'Interest Rate';'Load period';};
listUnits = {'';'kW';'$';'$/kW';'$/month';'%';'%';'years';};
listValues = cell(8,1);
listValues(1) = {testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Name};
listValues(2) = {num2str(testSystems(SYSINDEX).Generator(GENINDEX).Size)};
listValues(3) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost)};
listValues(4) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost/testSystems(SYSINDEX).Generator(GENINDEX).Size)};
listValues(5) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).OandM)};
listValues(6) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Financed)};
listValues(7) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanRate)};
listValues(8) = {num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanTerm)};
quote = '''';
for i=1:1:length(list)
    handles.(strcat('Cost_',list{i},'_text')) = uicontrol('Style', 'text', 'String', listText{i},...
    'Units','characters','Position', [50 49-5*i 35 1.5],'HorizontalAlignment','left',...
    'Tag', strcat('Cost_',list{i},'_text'),'FontSize', 12,'FontWeight','bold',...
    'Parent', handles.uipanelMain5','Visible','on');
    callback = strcat('@(hObject,eventdata)MainScreen1(',quote,'EditGenCost_Callback',quote,',hObject,eventdata,guidata(hObject))');
    handles.(strcat('Cost_',list{i},'_edit')) = uicontrol('Style', 'edit', 'String', listValues{i},...
    'Units','characters','Position', [50 47-5*i 20 2],'BackgroundColor',[1,1,1],...
    'Tag', strcat('Cost_',list{i},'_edit'),'FontSize', 10,...
    'Parent', handles.uipanelMain5','UserData',i,'Enable','on',...
    'Callback',eval(callback),'Visible','on');
    handles.(strcat('Cost_',list{i},'_unit')) = uicontrol('Style', 'text', 'String', listUnits{i},...
    'Units','characters','Position', [71 47.25-5*i 25 1.5],'HorizontalAlignment','left',...
    'Tag', strcat('Cost_',list{i},'_unit'),'FontSize', 10,...
    'Parent', handles.uipanelMain5','Visible','on');
end
handles.OptimizeGen = uicontrol('Style', 'pushbutton', 'String', 'Optimize Size',...
'Units','characters','Position', [80 36 30 3],'BackgroundColor',[1,0,0],...
'Tag', 'OptimizeGen','FontSize', 14,...
'Parent', handles.uipanelMain5','Callback','OptimizeGeneratorSize','Visible','on');

function EditGenCost_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
i = get(hObject,'UserData');
if i == 1
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Name = get(hObject,'String');
elseif i == 2
    testSystems(SYSINDEX).Generator(GENINDEX) = updateComponentSpec(testSystems(SYSINDEX).Generator(GENINDEX),'UB',str2double(get(hObject,'String')));
elseif i == 3
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost = str2double(get(hObject,'String'));
    set(handles.Cost_CostkW_edit,'String',num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost/testSystems(SYSINDEX).Generator(GENINDEX).Size));
elseif i == 4
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost = str2double(get(hObject,'String'))*testSystems(SYSINDEX).Generator(GENINDEX).Size;
    set(handles.Cost_Cost_edit,'String',num2str(testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Cost)); 
elseif i == 5
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).OandM = str2double(get(hObject,'String'));
elseif i == 6
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).Financed = str2double(get(hObject,'String'));
elseif i == 7
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanRate = str2double(get(hObject,'String'));
elseif i == 8
    testSystems(SYSINDEX).Costs.Equipment(GENINDEX).LoanTerm = str2double(get(hObject,'String'));
end
updateSummaryTable(handles)

function popupmenuBaseline_Callback(hObject, eventdata, handles)
%do nothing

% --- Executes during object creation, after setting all properties.
function popupmenuBaseline_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NPC_Years_Callback(hObject, eventdata, handles)
%Do nothing
function NPC_Years_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NPC_discount_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Costs.DiscountRate = str2double(get(hObject,'String'));

function NPC_discount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when entered data in editable cell(s) in uitableCosts.
function uitableCosts_CellEditCallback(hObject, eventdata, handles)
Data = get(hObject,'Data');
updateEqipmentCosts(Data);

function updateSummaryTable(handles)
global testSystems SYSINDEX
n = length(testSystems(SYSINDEX).Costs.Equipment);
Costs = cell(n,6);
for i = 1:1:n
    if ~isempty(testSystems(SYSINDEX).Costs.Equipment(i).Name)
        Costs(i,1) = {testSystems(SYSINDEX).Costs.Equipment(i).Name};
        Costs(i,2) = {testSystems(SYSINDEX).Costs.Equipment(i).Cost};
        Costs(i,3) = {testSystems(SYSINDEX).Costs.Equipment(i).OandM};
        Costs(i,4) = {testSystems(SYSINDEX).Costs.Equipment(i).Financed};
        Costs(i,5) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanRate};
        Costs(i,6) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm};
    end
end
set(handles.uitableCosts,'Data',Costs)

function updateEqipmentCosts(Data)
global testSystems SYSINDEX
for i = 1:1:length(Data(:,1))
    testSystems(SYSINDEX).Costs.Equipment(i).Name = Data{i,1};
    testSystems(SYSINDEX).Costs.Equipment(i).Cost = Data{i,2};
    testSystems(SYSINDEX).Costs.Equipment(i).OandM = Data{i,3};
    testSystems(SYSINDEX).Costs.Equipment(i).Financed = Data{i,4};
    testSystems(SYSINDEX).Costs.Equipment(i).LoanRate = Data{i,5};
    testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm = Data{i,6};
end

%% Tab 6 functions
function AggressiveOpt_Callback(hObject, eventdata, handles)
if get(handles.AggressiveOpt,'Value')
    set(handles.MedAncilOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function MedAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.MedAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.RobustAncilOpt,'Value',0)
end

function RobustAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.RobustAncilOpt,'Value')
    set(handles.AggressiveOpt,'Value',0)
    set(handles.MedAncilOpt,'Value',0)
end

function AutoAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.AutoAncilOpt,'Value')
    set(handles.ManualAncilOpt,'Value',0)
else
    set(handles.ManualAncilOpt,'Value',1)
end

function ManualAncilOpt_Callback(hObject, eventdata, handles)
if get(handles.ManualAncilOpt,'Value')
    set(handles.AutoAncilOpt,'Value',0)
else
    set(handles.AutoAncilOpt,'Value',1)
end

function uipanelOptimizationOptions_SelectionChangeFcn(hObject, eventdata, handles)
global testSystems SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'NoMixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = false;
    case 'MixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = true;
end

function Horizon_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Resolution_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));

function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on CarbonTax and none of its controls.
function CarbonTax_KeyPressFcn(hObject, eventdata, handles)


function DesignDay_Callback(hObject, eventdata, handles)


function StorageBuff_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
set(handles.editBuffer,'String',num2str(testSystems(SYSINDEX).Generator(stor(val)).VariableStruct.Buffer));

function StorageBuff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBuffer_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
val = get(handles.StorageBuff,'Value');
stor = get(handles.StorageBuff,'UserData');
testSystems(SYSINDEX).Generator(stor(val)).VariableStruct.Buffer = str2double(get(handles.editBuffer, 'String'));

function editBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function excessHeat_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessHeat = get(hObject, 'Value');

function excessCool_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessCool = get(hObject, 'Value');


% --- Executes on button press in SpinReserve.
function SpinReserve_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
if get(hObject,'Value')
    testSystems(SYSINDEX).optimoptions.SpinReserve = true;
    set(handles.SpinReservePerc,'Visible','on')
else
    testSystems(SYSINDEX).optimoptions.SpinReserve = false;
    set(handles.SpinReservePerc,'Visible','off')
end

function SpinReservePerc_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));

function SpinReservePerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
