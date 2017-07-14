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

% Last Modified by GUIDE v2.5 31-May-2017 23:42:54

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
Plant.optimoptions.forecast = 'Perfect';
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

files = dir(fullfile(Model_dir, 'Plant','*.mat'));
list=strrep({files.name},'.mat','');
set(handles.popupmenuProjectMain,'string',list)
set(handles.popupmenuProjectMain,'value',find(strcmp(testSystems(SYSINDEX).Name,list)))

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

if isfield(testSystems(SYSINDEX).Data,'Building')
    build = [testSystems(SYSINDEX).Data.Building.Type;testSystems(SYSINDEX).Data.Building.Climate;testSystems(SYSINDEX).Data.Building.Vintage;];
else build = 'custom';
end
set(handles.textBuilding,'String',build);
days = max(2,floor(testSystems(SYSINDEX).Data.Timestamp(end) - testSystems(SYSINDEX).Data.Timestamp(1)));
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
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/(days-1),10/(days-10)])
else
    set(handles.sliderDate1,'Min',0,'Max',days,'Value',0,'SliderStep',[1/(days-1),1/(days-1)])
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
            'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,'popupmenuAxes_Callback',quote,',hObject,eventdata,guidata(hObject))')),...
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
popupmenuAxes_Callback([], [], handles)


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
[f,p]=uiputfile(fullfile(Model_dir,'Plant','PlantNew.mat'),'Save Plant As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')

% --- Executes on selection change in popupmenuProjectMain.
function popupmenuProjectMain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant Model_dir testSystems SYSINDEX
% Load file that was selected from the popupmenu
SYSINDEX = length(testSystems)+1;
projList = get(handles.popupmenuProjectMain,'String');
projName = projList{get(handles.popupmenuProjectMain,'Value')};
projFile = fullfile(Model_dir,'Plant',projName);
load(projFile);
Plant.optimoptions.method = 'Planning';
Plant.optimoptions.forecast = 'Perfect';
allFieldNames = {'Name';'Data';'Generator';'optimoptions';'Network';'Costs';'subNet';'OpMatA';'OpMatB';'OneStep';'Online';'Design';'Dispatch';'Predicted';'RunData';'Baseline'};
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
if strcmp(n,'2')
    popupmenuDemand_Callback(hObject, eventdata, handles);
elseif strcmp(n,'1')
    set(handles.popupmenuAxes,'Value',1);
    popupmenuAxes_Callback(hObject, eventdata, handles)
end


%% Tab 1 functions

function System_Callback(hObject, eventdata, handles)
global SYSINDEX testSystems
% Change color
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
popupmenuAxes_Callback(hObject, eventdata, handles)
%change size/bold of correct button
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
Plant.optimoptions.forecast = 'Perfect';
%% go to system spec tab to edit


% --- Executes on button press in pushbuttonEDC.
function pushbuttonEDC_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Plant
Plant = testSystems(SYSINDEX);
% fNames = fieldnames(Plant);
% for i = 1:1:length(fNames)
%     if isempty(Plant.(fNames{i}))
%         Plant = rmfield(Plant,fNames{i});
%     end
% end
close
DISPATCH


% --- Executes on selection change in popupmenuAxes.
function popupmenuAxes_Callback(hObject, eventdata, handles)
global testSystems Plant RealTimeData TestData
list = get(handles.popupmenuAxes,'String');
item = list{get(handles.popupmenuAxes,'Value')};
if strcmp(item,'Monthly Costs')
    set(handles.sliderZoom1,'Visible','off');set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
    set(handles.textDay1,'Visible','off'); set(handles.textAllData1,'Visible','off'); set(handles.textHorizon1,'Visible','off');
    optim = get(get(handles.uipanelOptimizationOptions,'SelectedObject'),'Tag');
    handlesAll = guihandles;
    handlesAll.axesMain = handles.axesMain;
    handlesAll.axesCumulative = handles.axesCumulative;
    handles = handlesAll; %workaround to pass both axes handles and showSys handles
    for i_ts = 1:1:length(testSystems)
        if get(handles.(strcat('showSys',num2str(i_ts))),'Value') == 1
            if get(handles.DesignDay,'Value') ==1
                simTime = (testSystems(i_ts).Data.Timestamp(1):1/24*Plant.optimoptions.Resolution:testSystems(i_ts).Data.Timestamp(end))';
                [~, m, ~] = datevec(simTime);
                mUnique = unique(m); %months that be tested
                dStep = (testSystems(i_ts).Data.Timestamp(2) - testSystems(i_ts).Data.Timestamp(1));
                fullyear = min(365+1/24*testSystems(i_ts).optimoptions.Resolution,(testSystems(i_ts).Data.Timestamp(end)-testSystems(i_ts).Data.Timestamp(1) + dStep));
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
                if isempty(testSystems(i_ts).Design) || length(testSystems(i_ts).Design.Timestamp)~=fullyear*24/testSystems(i_ts).optimoptions.Resolution || any(testSystems(i_ts).Design.Timestamp==0) %at least some design days have not been run
                    repeat = true;
                else repeat = false;
                end
            end
            if  repeat
                Plant = testSystems(i_ts);
                TestData = Plant.Data;
                Xi = nnz(Plant.Data.Timestamp<=Plant.Data.Timestamp(1));
                Xf = nnz(Plant.Data.Timestamp<=Plant.Data.Timestamp(1)+365);
                RealTimeData = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);%create test data at correct frequency
                Plant.Design.Temperature = RealTimeData.Temperature;
                Plant.Design.Timestamp = zeros(length(RealTimeData.Temperature),1);
                Plant.Design.Demand = RealTimeData.Demand;
                if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
                    initializeOptimization%Load generators, build QP matrices
                end
                if ~isfield(Plant.Design,'GeneratorState') || length(Plant.Design.GeneratorState(:,1))~=length(Plant.Design.Timestamp)
                    % clear & initialize variables
                    nG = length(Plant.Generator);
                    [~,n] = size(Plant.OneStep.organize);
                    nL = n-nG;
                    Plant.Design.GeneratorState = zeros(length(Plant.Design.Timestamp),nG+nL);
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
    PlotHistorical(handles,item,1)
end

function RunPlanning(SiTest,interval)
global Plant Virtual RealTime DispatchWaitbar DateSim Last24hour 
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
        DateSim = Plant.Data.Timestamp(1) + (Si-1)*Plant.optimoptions.Resolution/24;
        TimeYesterday = DateSim-1+((0:round(24/Plant.optimoptions.Resolution)-1)'.*(Plant.optimoptions.Resolution/24));
        if Plant.Data.Timestamp(1)<(DateSim-1) && Plant.Data.Timestamp(end)>DateSim
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
        %change color and figure out spacing/width
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
DateSim = testSystems(SYSINDEX).Data.Timestamp(1) + get(handles.sliderDate1,'Value');
if floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))==0
    DateEnd = min(testSystems(SYSINDEX).Data.Timestamp(end),DateSim + 1);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))==1
    DateEnd = min(testSystems(SYSINDEX).Data.Timestamp(end),DateSim + 7);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))==2
    DateEnd = min(testSystems(SYSINDEX).Data.Timestamp(end),DateSim + 31);
elseif floor(get(handles.(strcat('sliderZoom',num2str(tab))),'Value'))==3
    DateEnd = min(testSystems(SYSINDEX).Data.Timestamp(end),DateSim + 365);
elseif get(handles.(strcat('sliderZoom',num2str(tab))),'Value')==4
    DateEnd = testSystems(SYSINDEX).Data.Timestamp(end);
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
popupmenuAxes_Callback(hObject, eventdata, handles)

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
popupmenuAxes_Callback(hObject, eventdata, handles)

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

% --- Executes on button press in Change.
function Change_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX Model_dir
%change building type
BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';};
Climate = {'Miami (ASHRAE 1A)';'Houston (ASHRAE 2A)';'Phoenix (ASHRAE 2B)';'Atlanta (ASHRAE 3A)';'Las Vegas (ASHRAE 3B-Inland)';'Los Angeles (ASHRAE 3B-Coast)';'San Francisco (ASHRAE 3C)';'Baltimore (ASHRAE 4A)';'Albuquerque (ASHRAE 4B)';'Seattle (ASHRAE 4C)';'Chicago (ASHRAE 5A)';'Boulder (ASHRAE 5B)';'Minneapolis (ASHRAE 6A)';'Helena, MT (ASHRAE 6B)';'Duluth, MN (ASHRAE 7)';'Fairbanks, AK (ASHRAE 8)';};
Vintage = {'2010 construction (ASHRAE 90.1-2010)';'2007 construction (ASHRAE 90.1-2007)';'2004 construction 90.1-2004';'“Post-1980” construction (ASHRAE 90.1-1989)';'“Pre-1980” construction';};
[s1,OK] = listdlg('PromptString','Select Building Type', 'SelectionMode','single','ListString',BuildType);
[s2,OK] = listdlg('PromptString','Select Building Type', 'SelectionMode','single','ListString',Climate);
[s3,OK] = listdlg('PromptString','Select Building Type', 'SelectionMode','single','ListString',Vintage);
testSystems(SYSINDEX).Data.Building.Type = BuildType{s1};
testSystems(SYSINDEX).Data.Building.Climate = Climate{s2};
testSystems(SYSINDEX).Data.Building.Vintage = Vintage{s3};
build = [testSystems(SYSINDEX).Data.Building.Type;testSystems(SYSINDEX).Data.Building.Climate;testSystems(SYSINDEX).Data.Building.Vintage;];
set(handles.textBuilding,'String',build)

vintageName = {'New2010';'New2007';'New2004';'Post1980';'Pre1980';};
climateName = {'_1A_'; '_2A_'; '_2B_'; '_3A_'; '_3B_'; '_3B-Coast_'; '_3C_'; '_4A_'; '_4B_'; '_4C_'; '_5A_'; '_5B_'; '_6A_'; '_6B_'; '_7_'; '_8_';};
buildTypeName = {'SDRest'; 'FFRest'; 'Sch-pri'; 'Sch-sec'; 'LgOff'; 'MdOff'; 'SmOff'; 'MRapt'; 'LgHotel'; 'SmHotel'; 'Hospital'; 'OutP'; 'Retail'; 'StMall'; 'SMarket'; 'ware';};
load(fullfile(Model_dir,'System Library','Buildings',strcat(buildTypeName(s1),climateName(s2),vintageName(s3)),'.mat'));
load(fullfile(Model_dir,'data','ClimateZoneTemperature.mat'));
%%% edit exterior lighting profile to avoid lowest points in energy profile
newMin (min(component.DemandE)+max(component.ExteriorLight));
AddLight = (component.ExteriorLight==0) & (component.DemandE<=newMin);
component.DemandE(AddLight) = newMin;
testSystems(SYSINDEX).Data.Timestamp = linspace(1/96,365,35040)'+datenum([2017,1,1]);
testSystems(SYSINDEX).Data.Temperature = Temperature(:,s2);
testSystems(SYSINDEX).Data.Demand.E = component.DemandE - component.CoolingElectricalLoad;%assume refrigeration loads are electric, but AC is under cooling demand
testSystems(SYSINDEX).Data.Demand.H = component.DemandH;
testSystems(SYSINDEX).Data.Demand.C = component.DemandC;


% --- Executes on selection change in popupmenuDemand.
function popupmenuDemand_Callback(hObject, eventdata, handles)
list = get(handles.popupmenuDemand,'String');
item = list{get(handles.popupmenuDemand,'Value')};
PlotHistorical(handles,item,2)

% --- Executes during object creation, after setting all properties.
function popupmenuDemand_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderDate2_Callback(hObject, eventdata, handles)
popupmenuDemand_Callback(hObject, eventdata, handles)

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
popupmenuDemand_Callback(hObject, eventdata, handles)

function sliderZoom2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderArea_Callback(hObject, eventdata, handles)

function sliderArea_CreateFcn(hObject, eventdata, handles)
function sliderWindowWall_Callback(hObject, eventdata, handles)

function sliderWindowWall_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderLighting_Callback(hObject, eventdata, handles)

function sliderLighting_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderPlugLoad_Callback(hObject, eventdata, handles)

function sliderPlugLoad_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderDensity_Callback(hObject, eventdata, handles)

function sliderDensity_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderUseHours_Callback(hObject, eventdata, handles)

function sliderUseHours_CreateFcn(hObject, eventdata, handles)
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

function compText4_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.eff = str2double(get(hObject,'String'));

function compText4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function compText3_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.Ramp = str2double(get(hObject,'String'));

function compText3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function compText2_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.LB = str2double(get(hObject,'String'));

function compText2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function compText1_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX GENINDEX
testSystems(SYSINDEX).Generator(GENINDEX).Size = str2double(get(hObject,'String'));
testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.UB = str2double(get(hObject,'String'));

function compText1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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

% --- Executes on selection change in CompFuel.
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
function popupmenuOptimization_Callback(hObject, eventdata, handles)

function popupmenuOptimization_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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


function StorageBuffer_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.Buffer = str2double(get(handles.editBuffer, 'String'));

function StorageBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function excessHeat_Callback(hObject, eventdata, handles)
global testSystems SYSINDEX
testSystems(SYSINDEX).optimoptions.excessHeat = get(hObject, 'Value');


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

%% Need to automatically create these control menus for all dispatchable generators
function popupmenuSys1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenuSys1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%Need to create extra sliders if there are more than 1 building
function sliderComfort_Callback(hObject, eventdata, handles)

function sliderComfort_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

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
        M(end+1) = d;
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


% --- Executes on key press with focus on CarbonTax and none of its controls.
function CarbonTax_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to CarbonTax (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

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
