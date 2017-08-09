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

% Last Modified by GUIDE v2.5 20-Jul-2017 08:04:33

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
global Plant Model_dir SYSINDEX GENINDEX testSystems 
movegui(gcf,'center');
set(gcf,'Name','Energy Planning Tool 2017.0.1')
Plant.optimoptions.method = 'Planning';
if isempty(Plant.Building)
    Plant.optimoptions.forecast = 'Perfect';
else Plant.optimoptions.forecast = 'Building';
end
if isempty(Plant.Costs) || ~isfield(Plant.Costs,'Equipment')
    defaultCosts
end
if isempty(testSystems)
    SYSINDEX = 1;%index in the list of systems (flull Plant structures)
    testSystems = Plant;
else
    Names = {};
    for i = 1:1:length(testSystems)
        Names(end+1) = {testSystems(i).Name};
    end
    SYSINDEX = find(strcmp(Plant.Name,Names));
    if isempty(SYSINDEX)
        SYSINDEX = length(testSystems)+1;
    end
    testSystems(SYSINDEX) = Plant;
end
GENINDEX = 1; %index in the list of generators of the current plant

%% Main Tabs
%Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2', etc.
TabText = {'Main Window';'Building Spec';'System Spec';'Network Spec';'Settings';};
for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        pos = get(handles.MainTab1','Position');
        set(handles.MainTab1,'Position',[pos(1),pos(2),pos(3),pos(4)+.5])
    else
        set(handles.(strcat('uipanelMain',j)),'Position',pan1pos)
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end
set(handles.axesMain,'NextPlot','add');
set(handles.axesCumulative,'NextPlot','add');
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

updateSystemRep(hObject, eventdata, handles)

files = dir(fullfile(Model_dir, 'Projects','*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuProjectMain,'string',list)
set(handles.popupmenuProjectMain,'value',find(strcmp(testSystems(SYSINDEX).Name,list)))

files = dir(fullfile(Model_dir, 'System Library','Buildings','*.mat'));
listB=strrep({files.name},'.mat','');
set(handles.popupmenuBuilding,'string',listB,'value',1);

files = dir(fullfile(Model_dir, 'System Library','Weather','*.mat'));
listW=strrep({files.name},'.mat','');
set(handles.popupmenuWeather,'string',listW,'value',1);
if ~isempty(testSystems(SYSINDEX).Building)
    I = find(strcmp(testSystems(SYSINDEX).Building.Name,listB));
    if isempty(I)
        listB(end+1) = {testSystems(SYSINDEX).Building.Name};
        I = length(listB);
    end
    set(handles.popupmenuBuilding,'string',listB);
    set(handles.popupmenuBuilding,'value',I);
    
    I = find(strcmp(testSystems(SYSINDEX).Weather.Name,listW));
    if isempty(I)
        listW(end+1) = {testSystems(SYSINDEX).Weather.Name};
        I = length(listW);
    end
    set(handles.popupmenuWeather,'string',listW);
    set(handles.popupmenuWeather,'value',I);
    
    set(handles.toggleHistorical,'Value',0);
    toggleSimulatedBuilding_Callback(hObject, eventdata, handles)
else
    toggleHistorical_Callback(hObject, eventdata, handles)
end

set(handles.sequential, 'value', testSystems(SYSINDEX).optimoptions.sequential);
set(handles.simultaneous, 'value', ~testSystems(SYSINDEX).optimoptions.sequential);
set(handles.excessHeat, 'value', testSystems(SYSINDEX).optimoptions.excessHeat);

set(handles.DesignDay,'value', 1);
set(handles.NoMixedInteger, 'value', ~testSystems(SYSINDEX).optimoptions.MixedInteger);
set(handles.MixedInteger, 'value', testSystems(SYSINDEX).optimoptions.MixedInteger);

set(handles.noSpinReserve, 'value', ~testSystems(SYSINDEX).optimoptions.SpinReserve);
set(handles.SpinReserve, 'value', testSystems(SYSINDEX).optimoptions.SpinReserve);
set(handles.SpinReservePerc, 'string', testSystems(SYSINDEX).optimoptions.SpinReservePerc);
set(handles.StorageBuffer, 'string', testSystems(SYSINDEX).optimoptions.Buffer);

n = length(testSystems(SYSINDEX).Costs.Equipment);
Costs = cell(n,6);
for i = 1:1:n
    Costs(i,1) = {testSystems(SYSINDEX).Costs.Equipment(i).Name};
    Costs(i,2) = {testSystems(SYSINDEX).Costs.Equipment(i).Cost};
    Costs(i,3) = {testSystems(SYSINDEX).Costs.Equipment(i).OandM};
    Costs(i,4) = {testSystems(SYSINDEX).Costs.Equipment(i).Financed};
    Costs(i,5) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanRate};
    Costs(i,6) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm};
end
set(handles.uitableCosts,'Data',Costs)
EditSystem(handles)
SysList_Make(handles)

if ~isempty(testSystems(SYSINDEX).Building)
    days = 365;
else
    days = max(2,floor(testSystems(SYSINDEX).Data.Timestamp(end) - testSystems(SYSINDEX).Data.Timestamp(1)));
end
if days<7
    set(handles.sliderZoom1,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data   
    set(handles.sliderZoom2,'Min',0,'Max',1,'Value',0,'SliderStep',[1,1]) %either single day or all data   
elseif days<31
    set(handles.sliderZoom1,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
    set(handles.sliderZoom2,'Min',0,'Max',2,'Value',0,'SliderStep',[1/2,1/2]) %either single day, week or all data
elseif days<367
    set(handles.sliderZoom1,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
    set(handles.sliderZoom2,'Min',0,'Max',3,'Value',0,'SliderStep',[1/3,1/3]) %either single day, week, month, or all data
else
    set(handles.sliderZoom1,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
    set(handles.sliderZoom2,'Min',0,'Max',4,'Value',0,'SliderStep',[1/4,1/4]) %either single day, week, month, year, or all data
end
if days>20
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
    set(handles.sliderDate2,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,10/days])
else
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
    set(handles.sliderDate2,'Min',0,'Max',days,'Value',0,'SliderStep',[1/days,1/days])
end
insertMockups(handles)

% --- Outputs from this function are returned to the command line.
function varargout = MainScreen1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function SysList_Make(handles)
global testSystems SYSINDEX
quote='''';
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
%Delete any old buttons
handlesAll = guihandles;
handlesAll.axesMain = handles.axesMain;
handlesAll.axesCumulative = handles.axesCumulative;
handles = handlesAll; %workaround to pass both axes handles and showSys handles
hNames = fieldnames(handles);
show = false(length(list),1);
for i = 1:1:length(hNames)
    if strncmp(hNames{i},'showSys',7) 
        if get(handles.(hNames{i}),'Value') == 1
            name = hNames{i};
            show(str2double(name(end))) = true;
        end
        delete(handles.(hNames{i}));
    end
    if strncmp(hNames{i},'System',6) || strncmp(hNames{i},'DeleteSys',9) || strncmp(hNames{i},'PrevSys',7) || strncmp(hNames{i},'NextSys',7) || strncmp(hNames{i},'textSys',7)
        delete(handles.(hNames{i}));
    end
end
show(SYSINDEX) = true;
%make new buttons
for tab = 1:1:4 %make gen list for tabs 1-4
    nSys = length(list);
    for i=1:1:nSys 
        num = num2str(i);
        curtab = floor((i-1)/3)+1;
        prev = 3*(curtab-1);
        xpos = 46 + 22*(i-prev);
        if tab == 3
            xpos = xpos + 40;
        end
        if curtab==1
            vis = 'on';
        else vis = 'off';
        end
        if i == SYSINDEX
            foregroundC = [0 0 0];
            backgroundC = [1 1 1];
            height = 2;
        else
            foregroundC = [0.5 0.5 0.5];
            backgroundC = [1 1 .7];
            height = 1.8;
        end
        handles.(strcat('System',num2str(tab),'_',num)) = uicontrol('Style', 'pushbutton', 'String', list{i},...
        'Units','characters',...
        'Position', [xpos 45 20 height],...
        'Tag', strcat('System',num2str(tab),'_',num),...
        'FontSize', 10,...
        'Parent', handles.(strcat('uipanelMain',num2str(tab))),...
        'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'System_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
        'Visible',vis,...
        'UserData',i,...
        'ForegroundColor', foregroundC,'BackgroundColor',backgroundC);
        if tab ==1 %Only make checkbox & delete options on Main Window
            %% checkbox to show
            handles.(strcat('showSys',num)) = uicontrol('Style', 'checkbox',...
            'Units','characters',...
            'Position', [xpos+10 43 3 1],...
            'Tag', strcat('showSys',num),...
            'Parent', handles.uipanelMain1,...
            'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'popupmenu_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
            'Visible',vis,'UserData',i);
            set(handles.(strcat('showSys',num)),'Value',show(i));
            %% text show:
            handles.(strcat('textSys',num)) = uicontrol('Style', 'text', 'String', 'Show:',...
            'Units','characters',...
            'Position', [xpos+1 43 8 1],...
            'Tag', strcat('textSys',num),...
            'Parent', handles.uipanelMain1,...
            'Visible',vis,'UserData',i);
            %% delete button
            handles.(strcat('DeleteSys',num)) = uicontrol('Style', 'pushbutton',...
            'Units','characters',...
            'Position', [xpos+15 42.5 5 2],...
            'Tag', strcat('DeleteSys',num),....
            'String','X','ForegroundColor','red','BackgroundColor',[.94,.94,.94],'FontSize',18,'FontWeight','bold',...
            'Parent', handles.uipanelMain1,...
            'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'DeleteSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
            'Visible',vis,'UserData',i);
        end
    end
    %make prev and next buttons
    handles.(strcat('PrevSys',num2str(tab))) = uicontrol('Style', 'pushbutton', 'String', 'Prev',...
    'Units','characters','FontSize', 10,'Position', [50 44 15 1.8],...
    'Tag', strcat('PrevSys',num2str(tab),'_',num),'Parent', handles.uipanelMain1,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'PrevSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
    'Visible','off','UserData',1);
    xpos = 140;
    if tab == 3
        xpos = 180;
    end
    handles.(strcat('NextSys',num2str(tab))) = uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Units','characters','FontSize', 10,'Position', [xpos 44 15 1.8],...
    'Tag', strcat('NextSys',num2str(tab),'_',num),'Parent', handles.uipanelMain1,...
    'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'NextSys_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
    'Visible','off','UserData',1);
    if nSys>3 
        set(handles.(strcat('NextSys',num2str(tab))),'Visible','on');
    end
end
popupmenu_Callback(handles.popupmenuAxes, [], handles)


% --- Executes on button press in PrevGen1.
function PrevSys_Callback(hObject, eventdata, handles)
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    tab=1;
elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
    tab=2;
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    tab=3;
elseif strcmp(get(handles.uipanelMain4,'Visible'),'on')
    tab=4;
end
list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.(strcat('PrevSys',num2str(tab))),'UserData');%current page of the list
if page<2
    set(handles.(strcat('PrevSys',num2str(tab))),'Visible','off','UserData',page-1)
else
    set(handles.(strcat('PrevSys',num2str(tab))),'Visible','on','UserData',page-1);
end
set(handles.(strcat('NextSys',num2str(tab))),'Visible','on','UserData',page-1)
for i = 1:1:12
    if 12*(page-1)+i<=gen
         j = num2str(12*(page-1)+i);
         set(handles.(strcat('System',num2str(tab),'_',j)),'Visible','off');
         if tab==1
             set(handles.(strcat('showSys',j)),'Visible','off');
             set(handles.(strcat('DeleteSys',j)),'Visible','off');
         end
    end
end
for i = 1:1:12
    j = num2str(12*(page-2)+i);
     set(handles.(strcat('System',num2str(tab),'_',j)),'Visible','on');
     if tab==1
         set(handles.(strcat('showSys',j)),'Visible','on');
         set(handles.(strcat('DeleteSys',j)),'Visible','on');
     end
end


% --- Executes on button press in NextGen1.
function NextSys_Callback(hObject, eventdata, handles)
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    tab=1;
elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
    tab=2;
elseif strcmp(get(handles.uipanelMain3,'Visible'),'on')
    tab=3;
elseif strcmp(get(handles.uipanelMain4,'Visible'),'on')
    tab=4;
end

list = get(handles.uipanelMain1,'UserData');
gen = length(list);
page = get(handles.(strcat('PrevSys',num2str(tab))),'UserData');%current page of the list
if page==ceil(gen/12)-1
    set(handles.(strcat('NextSys',num2str(tab))),'Visible','off','UserData',page+1)
else
    set(handles.(strcat('NextSys',num2str(tab))),'Visible','on','UserData',page+1);
end
set(handles.(strcat('PrevSys',num2str(tab))),'Visible','on','UserData',page+1)
for i = 1:1:12
     j = num2str(12*(page-1)+i);
     set(handles.(strcat('System',num2str(tab),'_',j)),'Visible','off');
     if tab==1
         set(handles.(strcat('showSys',j)),'Visible','off');
         set(handles.(strcat('DeleteSys',j)),'Visible','off');
     end
end
for i = 1:1:12
    j = num2str(12*(page)+i);
    if 12*(page)+i<=gen
         set(handles.(strcat('System',num2str(tab),'_',j)),'Visible','on');
         if tab==1
             set(handles.(strcat('showSys',j)),'Visible','on');
             set(handles.(strcat('DeleteSys',j)),'Visible','on');
         end
    end
end

% --- Executes on button press in saveProject.
function saveProject_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant Model_dir
Plant = testSystems(SYSINDEX);
[f,p]=uiputfile(fullfile(Model_dir,'Projects','New Project.mat'),'Save Project As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')

% --- Executes on selection savebuilding in popupmenuProjectMain.
function popupmenuProjectMain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant Model_dir testSystems SYSINDEX
% Load file that was selected from the popupmenu
SYSINDEX = length(testSystems)+1;
projList = get(handles.popupmenuProjectMain,'String');
projName = projList{get(handles.popupmenuProjectMain,'Value')};
projFile = fullfile(Model_dir,'Projects',projName);
load(projFile);
Plant.optimoptions.method = 'Planning';
if isempty(Plant.Building)
    Plant.optimoptions.forecast = 'Perfect';
else Plant.optimoptions.forecast = 'Building';
end
allFieldNames = {'Name';'Data';'Generator';'Building';'Weather';'Network';'Costs';'optimoptions';'subNet';'OpMatA';'OpMatB';'OneStep';'Online';'Design';'Dispatch';'Predicted';'RunData';'Baseline'};
fNames = fieldnames(Plant);
for i = 1:1:length(allFieldNames)
    if ~any(strcmp(allFieldNames{i},fNames))
        Plant.(allFieldNames{i}) = [];
    end
end
if isempty(Plant.Costs) || ~isfield(Plant.Costs,'Equipment')
    defaultCosts
end
testSystems(SYSINDEX) = Plant;
list={};
for i=1:length(testSystems)
    list(end+1) = {testSystems(i).Name};
end
set(handles.uipanelMain1,'UserData',list)
SysList_Make(handles);


% --- Executes during object creation, after setting all properties.
function popupmenuProjectMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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
if strcmp(n,'2')
    popupmenu_Callback(handles.popupmenuDemand, eventdata, handles);
elseif strcmp(n,'1')
    set(handles.popupmenuAxes,'Value',1);
    popupmenu_Callback(handles.popupmenuAxes, eventdata, handles)
end


%% Tab 1 functions

function System_Callback(hObject, eventdata, handles)
global SYSINDEX testSystems
% SaveBuilding color
for tab = 1:1:4
    pos = get(handles.(strcat('System',num2str(tab),'_',num2str(SYSINDEX))),'Position');
    pos(4) = 1.8;
    set(handles.(strcat('System',num2str(tab),'_',num2str(SYSINDEX))),'BackgroundColor',[1 1 .7],'ForegroundColor',[0.5 0.5 0.5],'Position',pos)
end
SYSINDEX = get(hObject,'UserData');
for tab = 1:1:4
    pos = get(handles.(strcat('System',num2str(tab),'_',num2str(SYSINDEX))),'Position');
    pos(4) = 2;
    set(handles.(strcat('System',num2str(tab),'_',num2str(SYSINDEX))),'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0],'Position',pos)
end
n = length(testSystems(SYSINDEX).Costs.Equipment);
Costs = cell(n,6);
for i = 1:1:n
    Costs(i,1) = {testSystems(SYSINDEX).Costs.Equipment(i).Name};
    Costs(i,2) = {testSystems(SYSINDEX).Costs.Equipment(i).Cost};
    Costs(i,3) = {testSystems(SYSINDEX).Costs.Equipment(i).OandM};
    Costs(i,4) = {testSystems(SYSINDEX).Costs.Equipment(i).Financed};
    Costs(i,5) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanRate};
    Costs(i,6) = {testSystems(SYSINDEX).Costs.Equipment(i).LoanTerm};
end
set(handles.uitableCosts,'Data',Costs)
popupmenu_Callback(hObject, eventdata, handles)
%% go to system spec tab to edit

function DeleteSys_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
SYSINDEX = get(hObject,'UserData');
if length(testSystems)>SYSINDEX
    testSystems(SYSINDEX:length(testSystems)-1) = testSystems(SYSINDEX+1:length(testSystems));
end
testSystems = testSystems(1:length(testSystems)-1);
SYSINDEX = 1;
SysList_Make(handles);


% --- Executes on button press in NewPlant.
function NewPlant_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.method = 'Planning';
if isempty(Plant.Building)
    Plant.optimoptions.forecast = 'Perfect';
else Plant.optimoptions.forecast = 'Building';
end
%% go to system spec tab to edit


% --- Executes on button press in pushbuttonEDC.
function pushbuttonEDC_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant
Plant = testSystems(SYSINDEX);
close
DISPATCH

function popupmenuAxes_Callback(hObject, eventdata, handles)
popupmenu_Callback(hObject, eventdata, handles)

function popupmenu_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant RealTimeData TestData
menu = char(get(hObject,'Tag'));
list = get(handles.(menu),'String');
item = list{get(handles.(menu),'Value')};
if strcmp(item,'Monthly Costs')
    set(handles.sliderZoom1,'Visible','off');set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
    set(handles.textDay1,'Visible','off'); set(handles.textAllData1,'Visible','off'); set(handles.textHorizon1,'Visible','off');
    handlesAll = guihandles;
    handlesAll.axesMain = handles.axesMain;
    handlesAll.axesCumulative = handles.axesCumulative;
    handles = handlesAll; %workaround to pass both axes handles and showSys handles
    for i_ts = 1:1:length(testSystems)
        if get(handles.(strcat('showSys',num2str(i_ts))),'Value') == 1
            if ~isempty(testSystems(i_ts).Building)
                fullyear = 365+1/24*testSystems(i_ts).optimoptions.Resolution;
            else
                dStep = (testSystems(i_ts).Data.Timestamp(2) - testSystems(i_ts).Data.Timestamp(1));
                fullyear = min(365+1/24*testSystems(i_ts).optimoptions.Resolution,(testSystems(i_ts).Data.Timestamp(end)-testSystems(i_ts).Data.Timestamp(1) + dStep));
            end
            if get(handles.DesignDay,'Value') ==1
                if ~isempty(testSystems(i_ts).Building)
                    mUnique = 1:1:12;
                else
                    simTime = (testSystems(i_ts).Data.Timestamp(1):1/24*testSystems(i_ts).optimoptions.Resolution:testSystems(i_ts).Data.Timestamp(end))';
                    [~, m, ~] = datevec(simTime);
                    mUnique = unique(m); %months that may need to be tested
                end
                SiTest = [];
                for j = 1:1:length(mUnique)
                    index = nonzeros((1:1:length(m))'.*(m==mUnique(j)));
                    SiTest(end+1) = index(1); %first point in each month
                end
                interval = 1;
                %% need to fix logical statements for when data is not 1 year
                if isempty(testSystems(i_ts).Design) || length(testSystems(i_ts).Design.Timestamp)~=(fullyear*24/testSystems(i_ts).optimoptions.Resolution) || any(testSystems(i_ts).Design.Timestamp(SiTest+1)==0) %at least some design days have not been run
                    repeat = true;
                else repeat = false;
                end
            else
                SiTest = 1;
                interval = fullyear;
                if isempty(testSystems(i_ts).Design) || length(testSystems(i_ts).Design.Timestamp)~=fullyear*24/testSystems(i_ts).optimoptions.Resolution || any(testSystems(i_ts).Design.Timestamp==0) %at least some points have not been run
                    repeat = true;
                else repeat = false;
                end
            end
            if  repeat
                Plant = testSystems(i_ts);
                TestData = [];
                if ~isempty(Plant.Data)
                    TestData = Plant.Data;
                    Xi = nnz(Plant.Data.Timestamp<=Plant.Data.Timestamp(1));
                    Xf = nnz(Plant.Data.Timestamp<=Plant.Data.Timestamp(1)+365);
                    RealTimeData = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);%create test data at correct frequency
                end
                if ~isfield(TestData,'Timestamp')
                    nS = 365*24/Plant.optimoptions.Resolution+1;
                    TestData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),nS)';
                end
                if ~isfield(TestData,'Temperature') && ~isempty(Plant.Weather)
                    D = datvec(TestData.Timestamp(1));
                    SoY = datenum([D(1),1,1]);
                    TestData.Temperature = interp1((0:1:8760)',[Plant.Weather.Tdb(1);Plant.Weather.Tdb],mod(24*(TestData.Timestamp-SoY),8760));
                end
                if ~isempty(Plant.Building)
                    %% add building data to extra columns in test data, or do nothing here?
                    if ~isfield(TestData,'Demand')
                        nS = length(TestData.Timestamp);
                        TestData.Demand.E = zeros(nS,1);
                        TestData.Demand.C = zeros(nS,1);
                        TestData.Demand.H = zeros(nS,1);
                    end
                    for i = 1:1:length(Plant.Building)
                        %use the network structure to find the location of the building
                        Location = [];
                        [Equipment,InteriorLighting,ExteriorLighting,Cooling,Heating,FanPower,OtherLoads] = BuildingProfile(Plant.Building(i),Plant.Weather,TestData.Timestamp,Location);
                        % Compare2Eplus(Plant.Building(i),Plant.Weather,TestData.Timestamp);
                        %% need conditions of if there are chillers and heaters only for electric
                        Electric = Equipment + InteriorLighting + ExteriorLighting + FanPower + Cooling/Plant.Building(i).VariableStruct.COP_C + Heating/Plant.Building(i).VariableStruct.COP_H + OtherLoads;
                        
                        TestData.Demand.E(:,end+1) = Electric;
                        TestData.Demand.C(:,end+1) = Cooling;
                        TestData.Demand.H(:,end+1) = Heating;
                        TestData.Holidays = [];
                    end
                    RealTimeData = TestData;
                end
                Plant.Design.Temperature = RealTimeData.Temperature;
                Plant.Design.Timestamp = zeros(length(RealTimeData.Temperature),1);
                Plant.Design.Demand = RealTimeData.Demand;
                if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
                    initializeOptimization%Load generators, build QP matrices
                end
                if ~isfield(Plant.Design,'GeneratorState') || length(Plant.Design.GeneratorState(:,1))~=length(Plant.Design.Timestamp)
                    Plant.Design.GeneratorState = zeros(length(Plant.Design.Timestamp),length(Plant.OneStep.organize));
                end
                RunPlanning(SiTest,interval);%specify starting indices Si, and duration of test in days
                testSystems(i_ts) = Plant;
                [testSystems(i_ts).Costs.Design,testSystems(i_ts).Costs.NPC]  = DesignCosts(i_ts);
            end
        end
    end
    %plot results
    PlotCosts(handles)
else
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
    if ~isempty(testSystems(SYSINDEX).Building)
        calculateBuildingProfile(handles,tab)
    else
        PlotHistorical(handles,item,tab)
    end
end

function RunPlanning(SiTest,interval)
global Plant Virtual RealTime DispatchWaitbar DateSim Last24hour TestData
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
Virtual = 1;
RealTime = 0;

Plant.optimoptions.Interval = interval;
NumOptim = ceil(Plant.optimoptions.Interval*24/Plant.optimoptions.Horizon);
timers = zeros(NumOptim*length(SiTest),3); % To record times set to zeros(1,3), to not record set to [];

if length(SiTest)==1
    STR = 'Optimizing Dispatch Throughout Entire Year';
else
    STR = 'Optimizing Design Day Dispatch';
end
DispatchWaitbar=waitbar(0,STR,'Visible','on');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
for j = 1:1:length(SiTest)
    Si=SiTest(j); %counter for # of times dispatch loop has run
    if isempty(Plant.Design) || Plant.Design.Timestamp(Si+1)==0
        DateSim = TestData.Timestamp(1) + (Si-1)*Plant.optimoptions.Resolution/24;
        TimeYesterday = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
        if TestData.Timestamp(1)<(DateSim-1) && TestData.Timestamp(end)>DateSim
            Last24hour = GetCurrentData(TimeYesterday);
        else %need to have this in terms of the first timestep
            Last24hour = GetCurrentData(TimeYesterday+1);
            Last24hour.Timestamp = TimeYesterday;
        end
        Plant.Design.GeneratorState(Si,:) = automaticInitialCondition(GetCurrentData(DateSim));
        LastDispatch = [];
        for t = 1:1:NumOptim
            Date = DateSim+[0;Time/24];
            Data = GetCurrentData(DateSim);
            Forecast = updateForecast(Date(2:end),Data);%% function that creates demand vector with time intervals coresponding to those selected
            [LastDispatch,timers((j-1)*NumOptim+t,:)] = DispatchLoop(Date,Forecast,LastDispatch);
        %     disp(strcat('FistDisp:',num2str(timers(Si,1))));
        %     disp(strcat('StebByStep:',num2str(timers(Si,2))));
        %     disp(strcat('FinalDisp:',num2str(timers(Si,3))));
            Si = StepDispatchForward(Si,Date,Data,Forecast,LastDispatch);
            waitbar(((j-1)*NumOptim+t)/(NumOptim*length(SiTest)),DispatchWaitbar,STR);
        end
        
    end
end
close(DispatchWaitbar)

function PlotCosts(handles)
global testSystems 
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
h1 = handles.axesMain;
cla(h1);
hold(h1,'on')
h2 = handles.axesCumulative;
cla(h2);
nShow = 0; %number of systems to plot
for i = 1:1:length(testSystems)
    if get(handles.(strcat('showSys',num2str(i))),'Value') == 1
        nShow = nShow+1;
    end
end
Xpos = zeros(12,nShow);
for i = 1:1:nShow
    Xpos(:,i) = (1/nShow^2 + (i-1)/nShow:1:12)';
end
j = 0;
NPC = zeros(nShow,1);
Name = cell(nShow,1);
colormap(h1,'summer')
for i = 1:1:length(testSystems)
    if get(handles.(strcat('showSys',num2str(i))),'Value') == 1
        j = j+1;
        Data = testSystems(i).Costs.Design;
        NPC(j) = testSystems(i).Costs.NPC;
        Name(j) = {testSystems(i).Name};
        %savebuilding color and figure out spacing/width
        bar(h1,Xpos(:,j),Data,'stacked','BarWidth',0.8/nShow)
    end
end
legend(h1,{'Financing Charges';'O & M Charges';'Demand Charges';'Electric Use Charges';'Fuel Charges';})
xlim(h1,[0,12])
ylabel(h1,'Cost ($)','Color','k','FontSize',14)
set(h1,'XTick',mean(Xpos,2),'XTickLabel', {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';})
colormap(h2,'summer')
bar(h2,NPC)
ylabel(h2,'20 Year Net Present Cost ($)','Color','k','FontSize',14)
set(h2,'XTickLabel',Name) 

function PlotHistorical(handles,demand,tab)
global testSystems SYSINDEX
if strcmp(get(handles.sliderZoom1,'Visible'),'off')
    set(handles.(strcat('sliderZoom',num2str(tab))),'Visible','on','value',1);set(handles.(strcat('sliderDate',num2str(tab))),'Visible','on','value',0);
    set(handles.(strcat('textDay',num2str(tab))),'Visible','on'); set(handles.(strcat('textAllData',num2str(tab))),'Visible','on'); set(handles.(strcat('textHorizon',num2str(tab))),'Visible','on');
end
%find the current date
DateSim = testSystems(SYSINDEX).Data.Timestamp(1) + get(handles.(strcat('sliderDate',num2str(tab))),'Value');
maxSlider = get(handles.(strcat('sliderZoom',num2str(tab))),'Max');
if isfield(testSystems(SYSINDEX).Data,'Demand')
    lastDate = testSystems(SYSINDEX).Data.Timestamp(end);
else lastDate = datenum([2018,1,1]);
end
if get(handles.(strcat('sliderZoom',num2str(tab))),'Value')==maxSlider
    DateEnd = lastDate;
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<1
    DateEnd = min(lastDate,DateSim + 1);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<2
    DateEnd = min(lastDate,DateSim + 7);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<3
    DateEnd = min(lastDate,DateSim + 31);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<4
    DateEnd = min(lastDate,DateSim + 365);
end
Xi = nnz(testSystems(SYSINDEX).Data.Timestamp<=DateSim);
Xf = nnz(testSystems(SYSINDEX).Data.Timestamp<=DateEnd);
time = testSystems(SYSINDEX).Data.Timestamp(Xi:Xf);
units = ' (kW)';
if strcmp(demand,'Electrical Demand')
    data = testSystems(SYSINDEX).Data.Demand.E(Xi:Xf,:);
    color = {'black'};
elseif strcmp(demand,'Cooling Demand')
    data = testSystems(SYSINDEX).Data.Demand.C(Xi:Xf,:);
    color = {'blue'};
elseif strcmp(demand,'Heating Demand')
    data = testSystems(SYSINDEX).Data.Demand.H(Xi:Xf,:);
    color = {'red'};
elseif strcmp(demand,'Water Demand')
    data = testSystems(SYSINDEX).Data.Demand.W(Xi:Xf,:);
    color = {'cyan'};
    units = ' (1000 CFS)';
end
if tab ==1 
    PlotData(handles.axesMain,time,data,[demand, units],color)
    PlotHistogram(handles.axesCumulative,data,[demand, units])
elseif tab ==2
    PlotData(handles.axesBuildLoad,time,data,[demand, units],color)
    PlotHistogram(handles.axesBuildCumulative,data,[demand, units])
end

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

function defaultCosts
global Plant
nG = length(Plant.Generator);
j = 0;
for i = 1:1:nG
    Type = Plant.Generator(i).Type;
    if isempty(strfind(Type,'Utility'))
        j = j+1;
        Plant.Costs.Equipment(j).Name = Plant.Generator(i).Name;
        if (strcmp(Type,'CHP Generator') || strcmp(Type,'Electric Generator')) && Plant.Generator(i).VariableStruct.isFuelCell
            costPerkW = 3000;
            OM = 100;
        elseif (strcmp(Type,'CHP Generator') || strcmp(Type,'Electric Generator')) 
            costPerkW = 1000;
            OM = 50;
        elseif strcmp(Type,'Chiller') %chillers
            if strcmp(Plant.Generator(i).Source,'Electricity')
                costPerkW = 100;
                OM = 10;
            else %absorption chiller
                costPerkW = 200;
                OM = 20;
            end
        elseif strcmp(Type,'Thermal Storage')%thermal storage
            costPerkW = 20;
            OM = 1;
        elseif strcmp(Type,'Electric Storage')%battery
            costPerkW = 500;
            OM = 10;
        elseif strcmp(Type,'Solar')%solar
            costPerkW = 500;
            OM = 5;
        end
        Plant.Costs.Equipment(j).Cost = costPerkW*Plant.Generator(i).Size;
        Plant.Costs.Equipment(j).OandM = OM*Plant.Generator(i).Size;
        Plant.Costs.Equipment(j).Financed = 100;
        Plant.Costs.Equipment(j).LoanRate = 6;
        Plant.Costs.Equipment(j).LoanTerm = 15;
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenuAxes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderDate1_Callback(hObject, eventdata, handles)
popupmenu_Callback(handles.popupmenuAxes, eventdata, handles)

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
popupmenu_Callback(handles.popupmenuAxes, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sliderZoom1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes when entered data in editable cell(s) in uitableCosts.
function uitableCosts_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableCosts (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
Data = get(hObject,'Data');
updateEqipmentCosts(Data);

%% Tab 2 functions

function popupmenuBuilding_Callback(hObject, eventdata, handles)
%load the selected building
global testSystems SYSINDEX Model_dir
list = get(handles.popupmenuBuilding,'string');
sel = get(handles.popupmenuBuilding,'Value');
load(fullfile(Model_dir,'System Library','Buildings',list{sel}));
testSystems(SYSINDEX).Building = building;
calculateBuildingProfile(handles,2)

function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuWeather_Callback(hObject, eventdata, handles)
%load selected weather profile
global testSystems SYSINDEX Model_dir
list = get(handles.popupmenuWeather,'string');
sel = get(handles.popupmenuWeather,'Value');
load(fullfile(Model_dir,'System Library','Weather',list{sel}));
testSystems(SYSINDEX).Weather = weather;
calculateBuildingProfile(handles,2)

function popupmenuWeather_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in toggleSimulatedBuilding.
function toggleSimulatedBuilding_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Model_dir
if isempty(testSystems(SYSINDEX).Building)
    list = get(handles.popupmenuWeather,'string');
    sel = get(handles.popupmenuWeather,'Value');
    load(fullfile(Model_dir,'System Library','Weather',list{sel}));
    testSystems(SYSINDEX).Weather = weather;
    popupmenuBuilding_Callback(hObject, eventdata, handles)
end
list = {'Electric';'Cooling';'Heating'};
set(handles.popupmenuDemand,'String',list);
list(end+1) = {'Monthly Costs'};
set(handles.popupmenuAxes,'String',list,'Value',1);
if get(handles.toggleHistorical,'Value')==1
    a = get(handles.toggleHistorical,'BackgroundColor');
    b = get(handles.toggleSimulatedBuilding,'BackgroundColor');
    c = get(handles.toggleHistorical,'ForegroundColor');
    d = get(handles.toggleSimulatedBuilding,'ForegroundColor');
    set(handles.toggleSimulatedBuilding,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.toggleHistorical,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.toggleSimulatedBuilding,'Value',1); %was already pressed
end
list = {'popupmenuBuilding';'popupmenuWeather';'SaveBuilding';'pushbuttonReset';'textDaily';'textInterHour';'textIntraHour';};
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','On');
end
set(handles.pushbuttonLoad,'Visible','Off')
pushbuttonPrev_Callback(hObject, eventdata, handles);
set(handles.editDailyVariability,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.DailyVariability))
set(handles.editHourlyVariability,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.HourlyVariability))
set(handles.editSubHourlyVariability,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.SubHourlyVariability))

% --- Executes on button press in toggleHistorical.
function toggleHistorical_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
if isempty(testSystems(SYSINDEX).Data)
    pushbuttonLoad_Callback(hObject, eventdata, handles)
end
if ~isfield(testSystems(SYSINDEX).Data,'Demand')
    uiwait(msgbox('Error: you must load demand data','Error','modal'))
else
    Outs =  fieldnames(testSystems(SYSINDEX).Data.Demand);
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
    else set(handles.toggleHistorical,'Value',1); %was already pressed
    end
    list = {'textBuilding';'pushbuttonPrev';'textArea';'textOccupancy';'textLighting';'textEquipment';'editComfort';'textComfort';'editRvalue';...
            'textClimate';'editArea';'editOccupancy';'pushbuttonOccupancy';'editLighting';'pushbuttonLighting';'editEquipment';'pushbuttonEquipment';...
            'pushbuttonNext';'textRvalue';'editAirChange';'textAirChange';'editDewPoint';'textDewPoint';'editColdSupply';...
            'textColdSupply';'editDamper';'textDamper';'popupmenuBuilding';'popupmenuWeather';'SaveBuilding';'pushbuttonReset';...
            'textVariability';'textDaily';'editDailyVariability';'textInterHour';'editHourlyVariability';'textIntraHour';'editSubHourlyVariability';};
    for i = 1:1:length(list)
        set(handles.(list{i}),'Visible','Off');
    end
    set(handles.pushbuttonLoad,'Visible','On')
end


function pushbuttonPrev_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = {'pushbuttonPrev';'editComfort';'textComfort';'editRvalue';'textRvalue';'editAirChange';'textAirChange';'editDewPoint';'textDewPoint';'editColdSupply';...
        'textColdSupply';'editDamper';'textDamper';'SaveBuilding';'pushbuttonReset';'editDailyVariability';'editHourlyVariability';'editSubHourlyVariability';};
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','Off');
end
set(handles.pushbuttonNext,'Visible','on');
set(handles.editArea,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.Area*10.76))
set(handles.textArea,'Visible','On')
set(handles.editOccupancy,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.occupancy*10.76))
set(handles.textOccupancy,'Visible','On')
set(handles.pushbuttonOccupancy,'Visible','On')
set(handles.editLighting,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.InteriorLights*10.76))
set(handles.textLighting,'Visible','On')
set(handles.pushbuttonLighting,'Visible','On')
set(handles.editEquipment,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.equipment*10.76))
set(handles.textEquipment,'Visible','On')
set(handles.pushbuttonEquipment,'Visible','On')

function pushbuttonNext_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
list = {'pushbuttonNext';'textArea';'textOccupancy';'textLighting';'textEquipment';'editArea';'editOccupancy';'pushbuttonOccupancy';'editLighting';'pushbuttonLighting';'editEquipment';'pushbuttonEquipment';};
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','Off');
end
set(handles.pushbuttonPrev,'Visible','On');
set(handles.editComfort,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.Comfort*9/5))
set(handles.textComfort,'Visible','On')
set(handles.editRvalue,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.Resistance*1055/3600*9/5*10.76/1000))
set(handles.textRvalue,'Visible','On')
set(handles.editAirChange,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.AirChangePerHr))
set(handles.textAirChange,'Visible','On')
set(handles.editDewPoint,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.DPset*9/5+32))
set(handles.textDewPoint,'Visible','On')
set(handles.editColdSupply,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.ColdAirSet*9/5+32))
set(handles.textColdSupply,'Visible','On')
set(handles.editDamper,'Visible','On','String',num2str(testSystems(SYSINDEX).Building.VariableStruct.MinDamper))
set(handles.textDamper,'Visible','On')


% --- Executes on button press in pushbuttonLoad.
function pushbuttonLoad_Callback(hObject, eventdata, handles)
%load historical data


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
calculateBuildingProfile(handles,2)

function editArea_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonOccupancy_Callback(hObject, eventdata, handles)
ScheduleEditor(hObject,'occ',handles)
calculateBuildingProfile(handles,2)

function editOccupancy_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.occupancy = str2double(get(hObject,'String'))/10.76;%convert to m^2
calculateBuildingProfile(handles,2)

function editOccupancy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonLighting_Callback(hObject, eventdata, handles)
ScheduleEditor(hObject,'light',handles)
calculateBuildingProfile(handles,2)

function editLighting_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.lighting = str2double(get(hObject,'String'))/10.76;%convert to m^2
calculateBuildingProfile(handles,2)

function editLighting_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonEquipment_Callback(hObject, eventdata, handles)
ScheduleEditor(hObject,'eqiuip',handles)
calculateBuildingProfile(handles,2)

function editEquipment_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.equipment = str2double(get(hObject,'String'))/10.76;%convert to m^2
calculateBuildingProfile(handles,2)

function editEquipment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editComfort_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.Comfort = str2double(get(hObject,'String'))*5/9;%convert to Celcius
calculateBuildingProfile(handles,2)

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
testSystems(SYSINDEX).Building.VariableStruct.Resistance = str2double(get(hObject,'String'))*3600/1055*5/9*1000/10.76;% convert hr-F-ft^2/BTU to m^2*K/kW
calculateBuildingProfile(handles,2)

function editRvalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAirChange_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.AirChangePerHr = str2double(get(hObject,'String'));
calculateBuildingProfile(handles,2)

function editAirChange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDewPoint_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.DPset = (str2double(get(hObject,'String'))-32)*5/9;
calculateBuildingProfile(handles,2)

function editDewPoint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editColdSupply_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.ColdAirSet = (str2double(get(hObject,'String'))-32)*5/9;
calculateBuildingProfile(handles,2)

function editColdSupply_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDamper_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.MinDamper = str2double(get(hObject,'String'));
calculateBuildingProfile(handles,2)

function editDamper_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDailyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.DailyVariability = str2double(get(hObject,'String'));
calculateBuildingProfile(handles,2)

function editDailyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editHourlyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.HourlyVariability = str2double(get(hObject,'String'));
calculateBuildingProfile(handles,2)

function editHourlyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSubHourlyVariability_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).Building.VariableStruct.SubHourlyVariability = str2double(get(hObject,'String'));
calculateBuildingProfile(handles,2)

function editSubHourlyVariability_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calculateBuildingProfile(handles,tab)
global testSystems SYSINDEX
%find the current date
DateSim = datenum([2017,1,1]) + get(handles.(strcat('sliderDate',num2str(tab))),'Value');
maxSlider = get(handles.(strcat('sliderZoom',num2str(tab))),'Max');
lastDate = datenum([2018,1,1]);
if get(handles.(strcat('sliderZoom',num2str(tab))),'Value')== maxSlider;
%     DateEnd = lastDate;
    A = datevec(DateSim);
    DateSim = datenum([A(1), 1, 1]);
    DateEnd = datenum([A(1)+1, 1, 1]);%1 year
    nS = 12*(DateEnd - DateSim); %2 hour INCREMENT
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<1
    DateEnd = min(lastDate,DateSim + 1);%1 day
    nS = 24*6; %10 MIN INCREMENT
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<2
    DateEnd = min(lastDate,DateSim + 7);%1 week
    nS = 24*7; %1 hour INCREMENT
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<3
    A = datevec(DateSim);
    DateSim = datenum([A(1), A(2), 1]);
    DateEnd = datenum([A(1), A(2)+1, 1]);%1 month
    nS = 24*(DateEnd - DateSim); %1 hour INCREMENT
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))<4
    A = datevec(DateSim);
    DateSim = datenum([A(1), A(2), 1]);
    DateEnd = datenum([A(1), min(12,A(2)+3), 1]);%3 month
    nS = 24*(DateEnd - DateSim); %1 hour INCREMENT
end
Date = linspace(DateSim,DateEnd,nS+1)';
nB = length(testSystems(SYSINDEX).Building);
Electric = zeros(nS+1,1);
for i = 1:1:nB
    [Equipment,InteriorLighting,ExteriorLighting,Cooling,Heating,FanPower,OtherLoads] = BuildingProfile(testSystems(SYSINDEX).Building(i),testSystems(SYSINDEX).Weather,Date,[]); % testSystems(SYSINDEX).Building(i).QPform.Location
    Electric = Electric + Equipment + InteriorLighting + ExteriorLighting + FanPower + Cooling/testSystems(SYSINDEX).Building(i).VariableStruct.COP_C + Heating/testSystems(SYSINDEX).Building(i).VariableStruct.COP_H + OtherLoads;
end
% Date = linspace(datenum([2017,1,1]),datenum([2018,1,1]),365*24/testSystems(SYSINDEX).optimoptions.Resolution+1)';
% Compare2Eplus(testSystems(SYSINDEX).Building(i),testSystems(SYSINDEX).Weather,Date,testSystems(SYSINDEX).Building(i).QPform.Location);
if tab ==1
    list = get(handles.popupmenuAxes,'String');
    item = list{get(handles.popupmenuAxes,'Value')};
elseif tab ==2
    list = get(handles.popupmenuDemand,'String');
    item = list{get(handles.popupmenuDemand,'Value')};
end
if strfind(item,'Electric')
    data = Electric;
    color = {'black'};
elseif strfind(item,'Cooling')
    data = Cooling;
    color = {'blue'};
elseif strfind(item,'Heating')
    data = Heating;
    color = {'red'};
end
units = ' (kW)';
if tab ==1
    PlotData(handles.axesMain,Date,data,[item, units],color)
    PlotHistogram(handles.axesCumulative,data,[item, units])
elseif tab ==2
    PlotData(handles.axesBuildLoad,Date,data,[item, units],color)
    PlotHistogram(handles.axesBuildCumulative,data,[item, units])
end

% --- Executes on selection savebuilding in popupmenuDemand.
function popupmenuDemand_Callback(hObject, eventdata, handles)
popupmenu_Callback(hObject, eventdata, handles)
% list = get(handles.popupmenuDemand,'String');
% item = list{get(handles.popupmenuDemand,'Value')};
% PlotHistorical(handles,item,2)

% --- Executes during object creation, after setting all properties.
function popupmenuDemand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderDate2_Callback(hObject, eventdata, handles)
popupmenu_Callback(handles.popupmenuDemand, eventdata, handles)

function sliderDate2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderZoom2_Callback(hObject, eventdata, handles)
if get(handles.sliderZoom2,'Value')== get(handles.sliderZoom2,'Max')
    set(handles.sliderDate2,'Visible','off');set(handles.textDate2,'Visible','off');
elseif strcmp(get(handles.sliderDate2,'Visible'),'off')
    set(handles.sliderDate2,'Visible','on');set(handles.textDate2,'Visible','on');
end
popupmenu_Callback(handles.popupmenuDemand, eventdata, handles)

function sliderZoom2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% Tab 3 functions
function SetupSystem(hObject,eventData,handles)
% This function updates the GUI graphic to represent the current plant
% configuration, showing buttons and lines connecting generating sources to
% the apropriate bus (AC, DC, thermal...)
global testSystems SYSINDEX GENINDEX
% Put Plant.Generator information into lists
nG = length(testSystems(SYSINDEX).Generator);
isType = zeros(nG,1);
isFC = zeros(nG,1);
type = cell(nG,1);
name = cell(nG,1);
source = cell(nG,1);
for i = 1:nG
    type(i) = {testSystems(SYSINDEX).Generator(i).Type};
    name(i) = {testSystems(SYSINDEX).Generator(i).Name};
    source(i) = {testSystems(SYSINDEX).Generator(i).Source};
    if (strcmp(type(i),'CHP Generator') || strcmp(type(i),'Electric Generator')) && testSystems(SYSINDEX).Generator(i).VariableStruct.isFuelCell
        isFC(i) = 1;
    end
end

% Find indices of relevant components
str = get(hObject,'String');
switch str
    case {'Fuel Cell';'ICE / mGT'}
        if strcmp(str,'Fuel Cell')
            isType = strcmp('CHP Generator',type).*isFC + strcmp('Electric Generator',type).*isFC;
        else
            isType = strcmp('CHP Generator',type).*(~isFC) + strcmp('Electric Generator',type).*(~isFC);
        end
    case 'Utility'
        isType = strcmp('Utility',type);
    case 'Solar PV'
        isType = strcmp('Solar',type);
    case 'Air Heater'
        isType = strcmp('Heater',type);
    case 'TES 2' % feeds Hot Water Demands
        isType = strcmp('Thermal Storage',type).*strcmp('Heat',source);
    case 'TES 3' % feeds Cooling Demands
        isType = strcmp('Thermal Storage',type).*strcmp('Cooling',source);
    case 'Battery'
        isType = strcmp('Electric Storage',type);
    case 'Heating Demands'
        GENINDEX = -1;
    case 'Hot Water Demands'
        GENINDEX = -2;
    case 'Cooling Demands'
        GENINDEX = -3;
    case 'AC / DC Conversion'
        GENINDEX = -4;
end
isType = nonzeros(linspace(1,nG,nG)'.*isType);

% Decide whether to show the user a list of the chosen type
if ~isempty(isType)
    if length(isType) ==1
        GENINDEX = isType;
    else
        % s - selection index; OK - whether an option was selected
        [s,OK] = listdlg('PromptString','Select desired component', ...
            'ListSize',[160 100], ...
            'SelectionMode','single', ...
            'ListString',name(isType));
        if OK
            GENINDEX = isType(s);
            EditSystem(handles)
        end
    end
end
EditSystem(handles)

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

%% The following need to represent different things for each system selected in 3rd tab
function CompName_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).Name = get(hObject,'String');


function CompName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function compText1_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(1,1) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'CHP Generator') || strcmp(Gen.Type,'Electric Generator')
    testSystems(SYSINDEX).Generator(GENINDEX).Size = str2double(get(hObject,'String'));
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.UB = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Size = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).Size = str2double(get(hObject,'String'));
end


function compText2_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(1,2) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'CHP Generator') || strcmp(Gen.Type,'Electric Generator')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.LB = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Sizem2 = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Voltage = str2double(get(hObject,'String'));
end

function compText3_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(2,1) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Eff = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MaxDOD = str2double(get(hObject,'String'));
end

function compText4_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(2,2) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Azimuth = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.ChargeResist = str2double(get(hObject,'String'));
end

function compText5_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(3,1) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Solar')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Tilt = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.DischResist = str2double(get(hObject,'String'));
end


function compText6_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Utility')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SumRates(3,2) = str2double(get(hObject,'String'));
elseif strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PeakCharge = str2double(get(hObject,'String'));
end

function compText7_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.PeakDisch = str2double(get(hObject,'String'));
end

function compText8_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Gen = testSystems(SYSINDEX).Generator(GENINDEX);
if strcmp(Gen.Type,'Electric Storage')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SelfDischarge = str2double(get(hObject,'String'))/(31*24*100);
end

function uitableEffCurve_CellEditCallback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
Outputs = fieldnames(testSystems(SYSINDEX).Generator(GENINDEX).Output);
nOutput = eventdata.Indices;
newValue = eventdata.NewData;
testSystems(SYSINDEX).Generator(GENINDEX).Output.(Outputs{nOutput(2)})(nOutput(1)) = newValue;

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
    emptyQPmatrices
end
updateSystemRep(hObject,[], handles)


% --- Executes on button press in pushbuttonRemove.
function pushbuttonRemove_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX 
if GENINDEX > 0
    str = strcat(testSystems(SYSINDEX).Generator(GENINDEX).Type,'.',testSystems(SYSINDEX).Generator(GENINDEX).Name);
    for n = 1:1:length(testSystems(SYSINDEX).Network)
        testSystems(SYSINDEX).Network(n).Equipment = testSystems(SYSINDEX).Network(n).Equipment(~strcmp(str,testSystems(SYSINDEX).Network(n).Equipment));
    end
    testSystems(SYSINDEX).Generator = testSystems(SYSINDEX).Generator([1:GENINDEX-1,GENINDEX+1:length(testSystems(SYSINDEX).Generator)]);
    GENINDEX = 0;
    updateSystemRep(hObject, eventdata, handles)
    EditSystem(handles)
    emptyQPmatrices
end

% --- Executes on selection savebuilding in CompFuel.
function CompFuel_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
fuels = get(hObject,'String');
testSystems(SYSINDEX).Generator(GENINDEX).Source = fuels{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function CompFuel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tab 5 functions
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
    case 'IdealDispatch'
        %Old DG-BEAT method
    case 'NoMixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = false;
        testSystems(SYSINDEX).optimoptions.method = 'Planning';
    case 'MixedInteger'
        testSystems(SYSINDEX).optimoptions.MixedInteger = true;
        testSystems(SYSINDEX).optimoptions.method = 'Planning';
    case 'NREL'
        testSystems(SYSINDEX).optimoptions.method = 'NREL';
end

function SpinningReserve_SelectionChangeFcn(hObject, eventdata, handles)
global testSystems SYSINDEX
switch get(eventdata.NewValue,'Tag')
    case 'noSpinReserve'
        testSystems(SYSINDEX).optimoptions.SpinReserve = false;
    case 'SpinReserve'
        testSystems(SYSINDEX).optimoptions.SpinReserve = true;
        testSystems(SYSINDEX).optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
end


function SpinReservePerc_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));

function SpinReservePerc_CreateFcn(hObject, eventdata, handles)
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


function Horizon_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function StorageBuffer_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Buffer = str2double(get(handles.editBuffer, 'String'));

function StorageBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on CarbonTax and none of its controls.
function CarbonTax_KeyPressFcn(hObject, eventdata, handles)



function excessHeat_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessHeat = get(hObject, 'Value');

function DesignDay_Callback(hObject, eventdata, handles)

function PlotHistogram(h,Data,Xlab)
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
cla(h)
bar(h,hist)
ylabel(h,'Percent of Time Within Range')
% set(h,'XTickLabel', label)
set(h,'XTickLabel','')
xlim(h,[0.5,length(hist)+.5]);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
pos = get(h,'position');
t = text(0,0,Xlab);
set(t,'Parent',h,'Units','characters','HorizontalAlignment','center','Position',[pos(3)/2,-6,0]);
% Place the text labels
for i = 1:length(hist)
    t = text(0,0,label{i});
    xpos = i*pos(3)/length(hist)- 0.5*pos(3)/length(hist);
    set(t,'Parent',h,'Units','characters','HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'Position',[xpos,0,0]);
end

function PlotData(h,TimeVec,Y,Ylab,color)
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
hold(h,'off')
for i = 1:1:n
    if i ==2
        hold(h,'on')
    end
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

function emptyQPmatrices
%something changed in the specifications, any previously calculated
%matrices or results do not apply
global testSystems SYSINDEX
A = {'subNet';'OpMatA';'OpMatB';'OneStep';'Online';'Dispatch';'Predicted';'RunData';'Baseline'};
for i = 1:1:length(A)
    testSystems(SYSINDEX).(A{i}) = [];
end

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

function checkboxSeasonal_Callback(hObject, eventdata, handles)
handlesc =get(handles.uipanelGenSpec,'Children');
if get(hObject,'Value')==1
    EditSystem('createPopup',handles,{'Summer','Winter'},[3 31.5 22 1],14,'bold','popupRates',1)
    x.Value = 1;
    MainScreen1('popupRates_Callback',x,0,handles)
    EditSystem('createText',handles,'From',[1 28.5 8 1.5],'Datetext',12,'bold')
    EditSystem('createText',handles,'To',[25 28.5 5 1.5],'Datetext',12,'bold')

else
    for i=1:length(handlesc)
        hide = {'popupRates','Datetext','popupDates1S','popupDates2S','popupDates3S','popupDates4S',...
            'popupDates1W','popupDates2W','popupDates3W','popupDates4W'};
        if nnz(strcmp(hide,get(handlesc(i),'Tag')))
            set(handlesc(i),'Visible','off')
            %delete(handlesc(i));
        end
    end
end

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


if ~isempty(strfind(Type,'S'))
    xS = xS(:,1:4);
    xS = fliplr(xS);
    if ~isempty(strfind(Type,'1'))
        Gen.VariableStruct.SumStartMonth = Date;
    elseif ~isempty(strfind(Type,'2'))
        Gen.VariableStruct.SumStartDay = Date;
    elseif ~isempty(strfind(Type,'3'))
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
    else
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
elseif ~isempty(strfind(Type,'W'))
    xW = xW(:,1:4);
    xW = fliplr(xW);
    if ~isempty(strfind(Type,'1'))
        Gen.VariableStruct.WinStartMonth = Date;
    elseif ~isempty(strfind(Type,'2'))
        Gen.VariableStruct.WinStartDay = Date;
    elseif ~isempty(strfind(Type,'3'))
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
    else
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

function createGridsellback(handles)
EditSystem('createText',handles,'Grid Sell Back',[2 25 30 1.5],'textEdit1',12,'bold')
EditSystem('createRadio',handles,'None',[10 23 15 1.5],10,'normal',1,'radiobuttonGridsellback')
EditSystem('createRadio',handles,'% of Tariffs',[10 21 19 1.5],10,'normal',0,'radiobuttonGridsellback')
EditSystem('createRadio',handles,'Reversed Meter',[10 18 22 1.5],10,'normal',0,'radiobuttonGridsellback')
EditSystem('createTextEdit',handles,'10',[14 19.5 8 1.5],'editTariffs',10,'normal','radiobuttonGridsellback')

              
function radiobuttonGridsellback_Callback(hObject,handles0,handles)
global GENINDEX SYSINDEX testSystems
handlesc = get(handles.uipanelGenSpec,'Children');
for i = 1:length(handlesc)
    if strcmp(get(handlesc(i),'Tag'),'radiobuttonGridsellback')
        if ~strcmp(get(handlesc(i),'String'),get(hObject,'String'))
            set(handlesc(i),'Value',0)
        else
            set(hObject,'Value',1)
        end
    elseif strcmp(get(handlesc(i),'Tag'),'editTariffs')
        x=i;
    end
end
if strcmp(get(hObject,'String'),'None')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = 0;
elseif strcmp(get(hObject,'String'),'% of Tariffs')
    editTariffs(hObject,handlesc(x))
elseif strcmp(get(hObject,'String'),'Reversed Meter')
    testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = -1;
end

function editTariffs(hObject,handles)
global GENINDEX SYSINDEX testSystems
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate = str2num(get(handles,'String'))/100;

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
elseif strcmp(get(hObject,'String'),'Bio-Fuel')
    testSystems(SYSINDEX).Generator(GENINDEX).Source = 'Biofuel';
else
    testSystems(SYSINDEX).Generator(GENINDEX).Source = get(hObject,'String');
end

function popupSolar_Callback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
stateNum = get(hObject,'Value');
state = stateName(stateNum);
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.State = char(state);

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

function uitableDCAC_CellEditCallback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Data = get(hObject,'Data');

function uitableBat_CellEditCallback(hObject,eventdata,handles)
global GENINDEX SYSINDEX testSystems
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.VoltCurve = get(hObject,'Data');


