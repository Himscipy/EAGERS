function EditSystem(varargin)
% This function updates the parameter options for the selected generator
% The user can then make changes to this particular generator and save
% those changes.
global testSystems SYSINDEX GENINDEX
if nargin && ischar(varargin{1})
    f= str2func(varargin{1});
    f(varargin{2:length(varargin)})
else
    editSystem(varargin{1})
end

function editSystem(handles)
global testSystems SYSINDEX GENINDEX
quote='''';
set(handles.uipanelLibrary,'Visible','off');
set(handles.uipanelGenSpec,'Visible','on');
set(handles.Library,'Visible','on');
if GENINDEX > 0
    Gen = testSystems(SYSINDEX).Generator(GENINDEX);
else
    switch GENINDEX
        case 0
            Gen = struct('Type', 'None', ...
                'Name', 'None');
        case -1
            Gen = struct('Type', 'Heating Demands', ...
                'Name', 'Heating Demands', ...
                'Demand', 100);
        case -2
            Gen = struct('Type', 'Hot Water Demands', ...
                'Name', 'Hot Water Demands', ...
                'Demand', 100);
        case -3
            Gen = struct('Type', 'Cooling Demands', ...
                'Name', 'Cooling Demands', ...
                'Demand', 100);
        case -4
            Gen = struct('Type', 'AC/DC Conversion', ...
                'Name', 'AC/DC Conversion', ...
                'Efficiency', 0.5);
    end
end
handlesC=get(handles.uipanelGenSpec,'Children');
for i= 1:length(handlesC)
    set(handlesC(i),'Visible','off')
    if strcmp(get(handlesC(i),'Tag'),'EffCurve')||strcmp(get(handlesC(i),'Tag'),'ResponseRate')
        clearGenAxes(handles)
    else
        delete(handlesC(i));
    end
end
createTextEdit(handles,Gen.Name,[28 34.5 75 2],'CompName',15,'bold')
switch Gen.Type
    case {'CHP Generator';'Electric Generator'}
        OutNames = {'Capacity';'Electricity';};
        Data(:,1) = Gen.Output.Capacity;
        Data(:,2) = Gen.Output.Electricity;
        editable = [true true];
        if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
            Data(:,3) = Gen.Output.Heat;
            OutNames = {'Capacity';'Electricity';'Heat';};
            editable = [true true true];
        end
        createTable(handles,Data,OutNames,'auto',{},[12 3 28.7 15.7143],'uitableEffCurve',editable)

        createText(handles,'Capacity (kW)',[70 32 20 1.75],'textEdit1',12,'normal')
        createText(handles,'Turn Down Ratio (x:1)',[60 30 30 1.75],'textEdit2',12,'normal') 
        createText(handles,'Response Rate',[6 32.5 22 1.75],'textEdit3',12,'normal')
        createText(handles,'Energy Source',[70 27 22 1.75],'textEdit4',12,'bold')
        
        createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal')
        createTextEdit(handles,num2str(Gen.Size/Gen.VariableStruct.Startup.Electricity(end)),[90 30 15 1.75],'compText2',10,'normal')
        v1=[];
        v2=[];
        v3=[];
        v4=[];
        v5=[];
        v6=[];
        if strcmp(Gen.Source, 'NG')
            v1=1;v2=0;v3=0;v4=0;v5=0;v6=0;    
        elseif strcmp(Gen.Source, 'Oil')
            v1=0;v2=1;v3=0;v4=0;v5=0;v6=0; 
        elseif strcmp(Gen.Source, 'Electricity')
            v1=0;v2=0;v3=1;v4=0;v5=0;v6=0; 
        elseif strcmp(Gen.Source, 'Heat')
            v1=0;v2=0;v3=0;v4=1;v5=0;v6=0; 
        elseif strcmp(Gen.Source, 'Steam')
            v1=0;v2=0;v3=0;v4=0;v5=1;v6=0; 
        elseif strcmp(Gen.Source, 'Biofuel')
            v1=0;v2=0;v3=0;v4=0;v5=0;v6=1; 
        end
                
        createRadio(handles,'Natural Gas',[65 25.5 22 1.75],10.5,'normal',v1,'radioSource')
        createRadio(handles,'Oil',[65 21.5 22 1.75],10.5,'normal',v2,'radioSource')
        createRadio(handles,'Electricity',[84 25.5 22 1.75],10.5,'normal',v3,'radioSource')
        createRadio(handles,'Heat',[84 21.5 22 1.75],10.5,'normal',v4,'radioSource')
        createRadio(handles,'Steam',[84 23.5 22 1.75],10.5,'normal',v5,'radioSource')
        createRadio(handles,'Bio-Fuel',[65 23.5 18 1.75],10.5,'normal',v6,'radioSource')
        
        %%Need to figure out plotting
        plotGenEfficiency(Gen,handles)
        plotResponse(Gen,handles)
    case 'Utility'
        if strcmp(Gen.Source,'Electricity')
            createText(handles,'(1) Off-Peak Rate ',[47 32 25 2],'textEdit1',12,'bold');
            createText(handles,'Electric Charge ($/kWh) ',[57 30 32 1.75],'textEdit2',11.5,'normal');
            createText(handles,'Demand Charge ($/kW) ',[57 28 32 1.75],'textEdit3',11.5,'normal');
            createText(handles,'(2) Partial-Peak Rate ',[47 26 31 1.75],'textEdit4',12,'bold');
            createText(handles,'Electric Charge ($/kWh) ',[57 24 32 1.75],'textEdit5',11.5,'normal');
            createText(handles,'Demand Charge ($/kW) ',[57 22 32 1.75],'textEdit6',11.5,'normal');
            createText(handles,'(3) Peak Rate ',[47 20 20 1.75],'textEdit7',12,'bold');
            createText(handles,'Electric Charge ($/kWh) ',[57 18 32 1.75],'textEdit8',11.5,'normal');
            createText(handles,'Demand Charge ($/kW) ',[57 16 32 1.75],'textEdit9',11.5,'normal');
            
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(1,1)),[90 30 15 1.75],'compText1',10,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(1,2)),[90 28 15 1.75],'compText2',10,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(2,1)),[90 24 15 1.75],'compText3',10,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(2,2)),[90 22 15 1.75],'compText4',10,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(3,1)),[90 18 15 1.75],'compText5',10,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.SumRates(3,2)),[90 16 15 1.75],'compText6',10,'normal');
            
            hrs = cell(24,1);
            for i = 1:1:length(hrs)
                hrs{i} = num2str(i);
            end
            colwidth = 21*ones(1,24);
            colwidth = num2cell(colwidth);
            
            createTable(handles,Gen.VariableStruct.SumRateTable,hrs,colwidth,{'Sun';'Mon';'Tue';'Wed';'Thu';'Fri';'Sat';},[2 1 102 10.5],'uitableEffCurve',true(1,24));

            if isfield(Gen.VariableStruct,'SumRateTable') && isfield(Gen.VariableStruct,'WinRateTable')
                handles = createCheckbox(handles,'Seasonal Rates',[1 33 18 1],1,'checkboxSeasonal');
                x.Value = 1;
                MainScreen1('checkboxSeasonal_Callback',handles.checkboxSeasonal,0,handles)
            else
                createCheckbox(handles,'Seasonal Rates',[1 33 18 1],0,'checkboxSeasonal');
            end
            MainScreen1('createGridsellback',handles)        
        else % Gen.Source: NG (Natural Gas)
            clearGenAxes(handles)
                               
            createText(handles,'Fuel Rate ($/MMBTU)',[60 32 30 1.75],'textEdit1',12,'normal');
            
            createTextEdit(handles,num2str(Gen.VariableStruct.Rate(1)),[90 32 15 1.75],'compText1',10,'normal');
        end
    case 'Heater'
        createText(handles,'Capacity (kW)',[70 32 20 1.75],'textEdit1',12,'normal');

        createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal');
    case 'Solar'
        createText(handles,'Location',[66 32 14 1.75],'textEdit1',12,'normal');
        createText(handles,'Size (kW)',[75 28 14 1.75],'textEdit2',12,'normal');
        createText(handles,'Size (m^2)',[73.5 26 16 1.75],'textEdit3',12,'normal');
        createText(handles,'Conversion Efficiency (%)',[55 22 34 1.75],'textEdit4',11.5,'normal');
        createText(handles,'Azimuth Angle (Degrees) ',[57 20 32 1.75],'textEdit5',11.5,'normal');
        createText(handles,'Tilt Angle (Degrees)',[63 18 26 1.75],'textEdit6',11.5,'normal');
        createText(handles,'Solar Type',[15 28.1 17 1.75],'textEdit7',12,'bold');
        createText(handles,'Solar Tracking',[12 24.1 22 1.75],'textEdit8',12,'bold');

        createTextEdit(handles,num2str(Gen.VariableStruct.Size),[90 28 15 1.75],'compText1',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Sizem2),[90 26 15 1.75],'compText2',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Eff),[90 22 15 1.75],'compText3',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Azimuth),[90 20 15 1.75],'compText4',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Tilt),[90 18 15 1.75],'compText5',10,'normal');

        stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';
         'Kentucky';'Louisiana';'Maine';'Maryland';'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'New Hampshire';'New Jersey';
         'New Mexico';'New York';'North Carolina';'North Dakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'Rhode Island';'South Carolina';'South Dakota';'Tennessee';'Texas';
         'Utah';'Vermont';'Virginia';'Washington';'West Virginia';'Wisconsin';'Wyoming';};
        stateNum = find(strcmp(Gen.VariableStruct.State,stateName));
        createPopup(handles,stateName,[80 32 25 1.75],10,'normal','popupSolar',stateNum);

        if strcmp(Gen.VariableStruct.PVtype,'flat')
            val1 = 1;
            val2 = 0;
        else
            val1 = 0;
            val2 = 1;
        end
        createRadio(handles,'Flat Panel',[5 27 16 1.5],11,'normal',val1,'radiobuttonSType');
        createRadio(handles,'Concentrated',[23 27 20 1.5],11,'normal',val2,'radiobuttonSType');

        if strcmp(Gen.VariableStruct.Tracking,'fixed')
            val1 = 1;
            val2 = 0;
            val3 = 0;
        elseif strcmp(Gen.VariableStruct.Tracking,'1axis')
            val1 = 0;
            val2 = 1;
            val3 = 0;
        else
            val1 = 0;
            val2 = 0;
            val3 = 1;
        end
        createRadio(handles,'Fixed',[1 23 10 1.5],11,'normal',val1,'radiobuttonSTracking');
        createRadio(handles,'Single Axis',[13 23 17 1.5],11,'normal',val2,'radiobuttonSTracking');            
        createRadio(handles,'Dual Axis',[32 23 15 1.5],11,'normal',val3,'radiobuttonSTracking');        

        row = {'DC rating';'Inverter/Transformer';'Mismatch';'Diodes/Connection';'DC wiring';...
               'AC wiring';'Soiling';'System availability';'Shading';'Sun-Tracking';'Age';};
        createTable(handles,Gen.VariableStruct.Data,{'Value';'Min';'Max'},'auto',row,[20 1.5 71.13 15.72],'uitableDCAC',logical([1 0 0]));

    case 'Electric Storage'
        createText(handles,'Size (kWh)',[70 32 20 1.75],'textEdit1',12,'normal');
        createText(handles,'Voltage (V)',[70 30 20 1.75],'textEdit2',12,'normal');
        createText(handles,'Max Depth of Discharge (%)',[69 27 20 3],'textEdit3',12,'normal');
        createText(handles,'Charging Internal Resistance (mOhms @ 100A)',[64 22 25 4],'textEdit4',12,'normal');
        createText(handles,'Discharging Internal Resistance  (mOhms @ 100A)',[61.5 16.5 28 4],'textEdit5',12,'normal');
        createText(handles,'Peak Charge Rate (C)',[9 32 31 1.75],'textEdit6',12,'normal');
        createText(handles,'Peak Disharge Rate (C)',[6 30 35 1.75],'textEdit7',12,'normal');
        createText(handles,'Self Discharge Rate (%)',[6.5 28 34 1.75],'textEdit8',12,'normal');

        createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Voltage),[90 30 15 1.75],'compText2',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.MaxDOD),[90 27.5 15 1.75],'compText3',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.ChargeResist),[90 23 15 1.75],'compText4',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.DischResist),[90 17.5 15 1.75],'compText5',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.PeakCharge),[41 32 15 1.75],'compText6',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.PeakDisch),[41 30 15 1.75],'compText7',10,'normal');
        createTextEdit(handles,num2str((Gen.VariableStruct.SelfDischarge)*(31*24*100)),[41 28 15 1.75],'compText8',10,'normal');

        createTable(handles,Gen.VariableStruct.VoltCurve,{'State of Charge';'Voltage'},'auto','numbered',[17 10 32.13 14.5],'uitableBat',true(1,2));

    case 'Thermal Storage'

    case 'Heating Demands'

    case 'Hot Water Demands'

    case 'Cooling Demands'

    case 'AC/DC Conversion'

    case 'None'
        
end
           
function plotGenEfficiency(Gen,handles)
axes(handles.EffCurve)
hold off
str = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity) > 0
    c = Gen.Output.Capacity./Gen.Output.Electricity;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Electricity, ...
        Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Electric'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
    if strcmp(Gen.Type,'CHP Generator')
        
       plot(AX(1),Gen.Output.Capacity,Gen.Output.Heat,'r-o');
    else
        c = Gen.Output.Capacity./Gen.Output.Heat;
        [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Heat, ...
            Gen.Output.Capacity(2:end),c(2:end));
        set(H1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
        hold on
    end
    str(end+1) = {'Heat'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam) > 0
    c = Gen.Output.Capacity./Gen.Output.Steam; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Steam, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','m','LineStyle','-','LineWidth',2,'Marker','o')
    hold on
    str(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling) > 0
    c = Gen.Output.Capacity./Gen.Output.Cooling; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Cooling, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','b','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Cooling'};
end
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','g','LineStyle',':','LineWidth',3)
str(end+1) = {'Cost Curve'};

if isfield(Gen.Output,'Cooling') && max(Gen.Output.Cooling) > 1
    ylim([0 max(Gen.Output.Cooling)])
    a = round(max(Gen.Output.Cooling))+1;
    ylim(AX(1),[0,a])
    ylim(AX(2),[0,1])
else
    a = round(max(c))+1;
    ylim(AX(1),[0,1])
    ylim(AX(2),[0,a])
    set(AX(1),'YTick',0:.1:1)
    set(AX(2),'YTick',0:a/10:a)
end
set(get(AX(1),'Ylabel'),'String','Efficiency')
set(get(AX(2),'Ylabel'),'String','Cost Curve Shape')
xlabel('% of Capacity')
legend(str);
title('Efficiency / Cost')
set(handles.EffCurve,'UserData',AX)
set(handles.EffCurve,'Tag','EffCurve')
set(AX(2),'Tag','EffCurve')


function plotResponse(Gen,handles)
A = Gen.VariableStruct.StateSpace.A;
B = Gen.VariableStruct.StateSpace.B;
C = Gen.VariableStruct.StateSpace.C;
D = Gen.VariableStruct.StateSpace.D;
Dt = Gen.VariableStruct.StateSpace.Dt;
SS = ss(A,B,C,D,Dt);
x0 = [];
Names = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Electricity(end);
    Names(end+1) = {'Electricity'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling) > 0

    x0(end+1) = Gen.VariableStruct.Startup.Cooling(end);
    Names(end+1) = {'Cooling'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Steam(end);
    Names(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Heat(end);
    Names(end+1) = {'Heat'};
end
nS = round(3600/Dt)+1;
t = linspace(0, Dt*(nS-1),nS);
u = Gen.Size*linspace(1,1,nS);
[n,n2] = size(C);
X0 = zeros(n2,1);
for i = 1:1:n
    X0(find(C(i,:),1,'first'))=x0(i);
end
% r = size(SS(1,1).A);
% s = length(x0);
% x0 = [x0;zeros(r(1)-s,1);];
[y,t] = lsim(SS,u,t,X0);
RR = plot(handles.ResponseRate,t/60,y);
xlabel('Time (min)')
legend(Names)

set(handles.ResponseRate,'UserData',RR)
set(handles.ResponseRate,'Tag','ResponseRate')

function clearGenAxes(handles)
AX_eff = get(handles.EffCurve,'UserData');
AX_res = get(handles.ResponseRate,'UserData');
if ~isempty(AX_eff)
    cla(AX_eff(1))
    cla(AX_eff(2))
    set(AX_eff(2),'Visible','off')
    legend(AX_eff(1),'hide')
end
if ~isempty(AX_res)
    set(AX_res,'Visible','off')
end

function createTable(handles,Data,ColName,ColWidth,RowName,position,tag,edit)
quote='''';
Table = uitable('Parent',handles.uipanelGenSpec,...
                 'Units','characters',...
                 'Position',position,...
                 'FontSize',8,...
                 'Data',Data,...
                 'RowName',RowName,...
                 'ColumnName',ColName,...
                 'ColumnWidth',ColWidth,...
                 'Tag',tag,...
                 'ColumnEditable',edit,...
                 'CellEditCallback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,tag,'_CellEditCallback',quote,',hObject,eventdata,guidata(hObject))')));
              
function createText(handles,string,position,tag,size,weight)
Text = uicontrol('Style','text',...
                 'Parent', handles.uipanelGenSpec,...
                 'Units','characters',...
                 'Position',position,...
                 'String',string,...
                 'FontSize',size,...
                 'FontWeight',weight,...
                 'Tag',tag);
             
function createTextEdit(handles,string,position,tag,size,weight,varargin)
quote='''';
if ~isempty(varargin)%if the tag and call are different
    call = varargin{1};
else
    call = tag;
end
textEdit = uicontrol('Style','edit',...
                 'Parent',handles.uipanelGenSpec,...
                 'String',string,...
                 'Units','characters',...
                 'Position',position,...
                 'FontSize',size,...
                 'FontWeight',weight,...
                 'Tag',tag,...
                 'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,call,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

function handles = createCheckbox(handles,string,position,value,tag)
quote='''';
handles.(tag) = uicontrol('Style','checkbox','String',string,...
                     'Parent',handles.uipanelGenSpec,...
                     'Units','Characters',...
                     'Position',position,...
                     'Value',value,...
                     'Tag',tag,...
                     'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,tag,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

function createPopup(handles,string,position,size,weight,tag,value,varargin)
quote='''';
if ~isempty(varargin)%if the tag and call are different
    call = varargin{1};
else
    call = tag;
end
Popup = uicontrol('Style','popup',...
                       'String',string,...
                       'Value',value,...
                       'Parent',handles.uipanelGenSpec,...
                       'Units','Characters',...
                       'Position',position,...
                       'Tag',tag,...
                       'FontSize',size,...
                       'FontWeight',weight,...
                       'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,call,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

function createRadio(handles,string,position,size,weight,value,tag)
quote='''';
Radio = uicontrol('Style','radiobutton',...
                  'String',string,...
                  'Parent',handles.uipanelGenSpec,...
                  'Units','Characters',...
                  'Position',position,...
                  'Value',value,...
                  'FontSize',size,...
                  'FontWeight',weight,...
                  'Tag',tag,...
                  'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,tag,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

