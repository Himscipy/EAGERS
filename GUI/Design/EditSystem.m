function EditSystem(varargin)
% This function updates the parameter options for the selected generator
% The user can then make changes to this particular generator and save
% those changes.
if nargin && ischar(varargin{1})
    f= str2func(varargin{1});
    f(varargin{2:length(varargin)})
else
    editSystem(varargin{1})
end

function editSystem(handles)
global testSystems SYSINDEX GENINDEX
set(handles.uipanelLibrary,'Visible','off');
set(handles.uipanelGenSpec,'Visible','on');
set(handles.Library,'Visible','on');
set(handles.saveSystem,'Visible','on');
set(handles.pushbuttonRemove,'Visible','on');
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
    delete(handlesC(i));
end
createTextEdit(handles,Gen.Name,[28 34.5 75 2],'CompName',15,'bold')
switch Gen.Type
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
            handles = createCheckbox(handles,'Seasonal Rates',[1 33 18 1],1,'checkboxSeasonal');
            EditSystem('checkboxSeasonal_Callback',handles.checkboxSeasonal,[],handles)
            EditSystem('createText',handles,'Grid Sell Back',[2 24.5 30 1.5],'textEdit10',12,'bold')
            
            if testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate == -1
                v1=0;v2=1;v3=0;
            elseif testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh>=0 && testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate==0
                v1=1;v2=0;v3=0;
            else
                v1=0;v2=0;v3=1;
            end
        
            EditSystem('createRadio',handles,'None',[2 23 10 1.5],10,'normal',v1,'None','radiobuttonGridsellback')
            EditSystem('createRadio',handles,'% Purchase',[2 21 18 1.5],10,'normal',v2,'SellbackAsPerc','radiobuttonGridsellback')
            EditSystem('createRadio',handles,'Fixed Rate',[2 19 18 1.5],10,'normal',v3,'FixedRate','radiobuttonGridsellback')
            if v1 == 1
                createText(handles,'Minimum Import (kW)',[10 16 30 1.75],'textEdit11',12,'normal');
            else
                createText(handles,'Maximum Export (kW)',[10 16 30 1.75],'textEdit11',12,'normal');
            end
            EditSystem('createTextEdit',handles,num2str(abs(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.MinImportThresh)),[20 14 12 1.5],'maxSellback',10,'normal')
            EditSystem('createTextEdit',handles,num2str(testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackPerc),[20 21 8 1.5],'editTariffs',10,'normal')
            EditSystem('createTextEdit',handles,num2str(max(0,testSystems(SYSINDEX).Generator(GENINDEX).VariableStruct.SellBackRate)),[20 19 8 1.5],'editSellbackRate',10,'normal')
            createText(handles,'$/kWh',[28 18.75 10 1.75],'textEdit12',11.5,'normal');
        else % Gen.Source: NG (Natural Gas)
            createText(handles,'Fuel Rate ($/MMBTU)',[60 32 30 1.75],'textEdit1',12,'normal');
            createTextEdit(handles,num2str(Gen.VariableStruct.Rate(1)),[90 32 15 1.75],'compText1',10,'normal');
        end
    case {'CHP Generator';'Electric Generator';'Heater';'Chiller';}
        if strcmp(Gen.Type,'Heater')
            OutNames = {'Capacity';'Heat';};
            Data = [Gen.Output.Capacity,Gen.Output.Heat;];
            pos = [16 3 19.3 15.7143];
            LB = Gen.VariableStruct.Startup.Heat(end);
        elseif strcmp(Gen.Type,'Chiller')
            OutNames = {'Capacity';'Cooling';};
            pos = [10 3 20 16];
            Data = [Gen.Output.Capacity,Gen.Output.Cooling;];
            LB = Gen.VariableStruct.Startup.Cooling(end);
        elseif strcmp(Gen.Type,'Electric Generator')
            OutNames = {'Capacity';'Electricity';};
            pos = [10 3 20 16];
            Data = [Gen.Output.Capacity,Gen.Output.Electricity;];
            LB = Gen.VariableStruct.Startup.Electricity(end);
        elseif strcmp(Gen.Type,'CHP Generator') 
            pos = [8 3 30 16];
            Data = [Gen.Output.Capacity,Gen.Output.Electricity,Gen.Output.Heat;];
            OutNames = {'Capacity';'Electricity';'Heat';};
            LB = Gen.VariableStruct.Startup.Electricity(end);
        end
        createTable(handles,Data,OutNames,'auto',{},pos,'uitableEffCurve',true(1,length(Data(1,:))))
        
        createText(handles,'Capacity (kW)',[70 32 20 1.75],'textEdit1',12,'normal')
        createText(handles,'Minimum Output (kW)',[60 30 30 1.75],'textEdit2',12,'normal')
        createText(handles,'Ramp Rate (kW / hr)',[60 28 30 1.75],'textEdit3',12,'normal') 
        createText(handles,'Startup Cost ($/start)',[60 26 30 1.75],'textEdit4',12,'normal') 
        createText(handles,'Nat. Freq.',[1 32.5 15 1.65],'textEdit5',12,'normal')  
        createText(handles,'Damping',[32 32.5 13 1.65],'textEdit6',12,'normal')  
        createText(handles,'Energy Source',[70 23.5 22 1.75],'textEdit7',12,'bold')
        
        createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal')
        createTextEdit(handles,num2str(LB),[90 30 15 1.75],'compText2',10,'normal')
        if isfield(Gen.VariableStruct,'dX_dt')
            createTextEdit(handles,num2str(Gen.VariableStruct.dX_dt),[90 28 15 1.75],'compText3',10,'normal')
            p = eig(Gen.VariableStruct.StateSpace.A);
            w_0 = sqrt(real(p(1))^2 + imag(p(1))^2);
            zeta = -real(p(1)+p(2))/(2*w_0);
            createTextEdit(handles,num2str(w_0),[16 32.5 14 1.75],'compText5',10,'normal')
            createTextEdit(handles,num2str(zeta),[45 32.5 14 1.75],'compText6',10,'normal')
        else
            createTextEdit(handles,'--',[90 28 15 1.75],'compText3',10,'normal')
            createTextEdit(handles,'--',[16 32.5 14 1.75],'compText5',10,'normal')
            createTextEdit(handles,'--',[45 32.5 14 1.75],'compText6',10,'normal')
        end
        if isfield(Gen.VariableStruct,'StartCost')
            createTextEdit(handles,num2str(Gen.VariableStruct.StartCost),[90 26 15 1.75],'compText4',10,'normal')
        else
            createTextEdit(handles,'--',[90 26 15 1.75],'compText4',10,'normal')
        end
        
        if strcmp(Gen.Source, 'NG')
            v1=1;v2=0;v3=0;v4=0;
        elseif strcmp(Gen.Source, 'Oil')
            v1=0;v2=1;v3=0;v4=0;
        elseif strcmp(Gen.Source, 'Electricity')
            v1=0;v2=0;v3=1;v4=0;
        elseif strcmp(Gen.Source, 'Heat')
            v1=0;v2=0;v3=0;v4=1;
        end
        
        if strcmp(Gen.Type,'Heater')
            createRadio(handles,'Natural Gas',[60 22 25 1.75],10.5,'normal',v1,'NatGas','radioSource')
%             createRadio(handles,'Oil',[87 22 18 1.75],10.5,'normal',v2,'radioSource')
            createRadio(handles,'Electricity',[85 22 18 1.75],10.5,'normal',v3,'Electricity','radioSource')
        elseif strcmp(Gen.Type,'Chiller')
            createRadio(handles,'Electricity',[60 22 25 1.75],10.5,'normal',v3,'Electricity','radioSource')
            createRadio(handles,'Heat',[85 22 18 1.75],10.5,'normal',v4,'Heat','radioSource')
        else
            createRadio(handles,'Natural Gas',[60 22 25 1.75],10.5,'normal',v1,'NatGas','radioSource')
            createRadio(handles,'Diesel',[85 22 18 1.75],10.5,'normal',v2,'Diesel','radioSource')
        end
            
        %%Need to figure out plotting
        
        handles.ResponseRate = axes('Units','characters',...
            'Position', [9,22,40,10],'NextPlot','add',...
            'Tag', 'ResponseRate',...
            'Parent', handles.uipanelGenSpec,...
            'Visible','on');

        handles.EffCurve = axes('Units','characters',...
            'Position', [58,3.5,44,15],'NextPlot','add',...
            'Tag', 'EffCurve',...
            'Parent', handles.uipanelGenSpec,...
            'Visible','on');

        plotGenEfficiency(Gen,handles)
        secondOrderResponse(Gen,handles);

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
        createRadio(handles,'Flat Panel',[5 27 16 1.5],11,'normal',val1,'Flat','radiobuttonSType');
        createRadio(handles,'Concentrated',[23 27 20 1.5],11,'normal',val2,'Concentrated','radiobuttonSType');

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
        createRadio(handles,'Fixed',[1 23 10 1.5],11,'normal',val1,'Fixed','radiobuttonSTracking');
        createRadio(handles,'Single Axis',[13 23 17 1.5],11,'normal',val2,'SingleAxis','radiobuttonSTracking');            
        createRadio(handles,'Dual Axis',[32 23 15 1.5],11,'normal',val3,'DualAxis','radiobuttonSTracking');        

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
        createText(handles,'Size (kWh)',[70 32 20 1.75],'textEdit1',12,'normal');
        createText(handles,'Size (L)',[72 30 20 1.75],'textEdit2',12,'normal');
        createText(handles,'T Hot(C)',[72 28 20 1.75],'textEdit3',12,'normal');
        createText(handles,'T Cold(C)',[71.5 26 20 1.75],'textEdit4',12,'normal');
        createText(handles,'Ramp Rate',[71 24 20 1.75],'textEdit5',12,'normal');
        createText(handles,'Charging Efficiency (%)',[10.75 32 32 1.75],'textEdit6',12,'normal');
        createText(handles,'Discharging Efficiency (%)',[6 30 38 1.75],'textEdit7',12,'normal');
        createText(handles,'Self Discharge(% lost per day)',[1.5 28 42 1.75],'textEdit8',12,'normal');
        createText(handles,'Fill Rate (L/min)',[15 26 34 1.75],'textEdit9',12,'normal');
        createText(handles,'Discharge Rate (L/min)',[10.25 24 34 1.75],'textEdit10',12,'normal');
        
        createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.SizeLiter),[90 30 15 1.75],'compText2',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Tcold),[90 28 15 1.75],'compText3',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.Thot),[90 26 15 1.75],'compText4',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.dX_dt),[90 24 15 1.75],'compText5',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.ChargeEff),[45 32 15 1.75],'compText6',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.DischargeEff),[45 30 15 1.75],'compText7',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.SelfDischarge),[45 28 15 1.75],'compText8',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.FillRate),[45 26 15 1.75],'compText9',10,'normal');
        createTextEdit(handles,num2str(Gen.VariableStruct.DischRate),[45 24 15 1.75],'compText10',10,'normal');
        
%     case 'Boiler'
%         OutNames = {'Capacity';'Steam';};
%         Data(:,1) = Gen.Output.Capacity;
%         Data(:,2) = Gen.Output.Steam;
%         editable = [true true];
%         createTable(handles,Data,OutNames,'auto',{},[16 3 19.3 15.7143],'uitableEffCurve',editable)
%         
%         createText(handles,'Capacity (kW)',[70 32 20 1.75],'textEdit1',12,'normal')
%         
%         createTextEdit(handles,num2str(Gen.Size),[90 32 15 1.75],'compText1',10,'normal')
%         
%         plotGenEfficiency(Gen,handles)
        
    case 'Heating Demands'

    case 'Hot Water Demands'

    case 'Cooling Demands'

    case 'AC/DC Conversion'

    case 'None'
        
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
                     'Callback',eval(strcat('@(hObject,eventdata)EditSystem(',quote,tag,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

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

function createRadio(handles,string,position,size,weight,value,tag,call)
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
                  'Callback',eval(strcat('@(hObject,eventdata)MainScreen1(',quote,call,'_Callback',quote,',hObject,eventdata,guidata(hObject))')));

function checkboxSeasonal_Callback(hObject,eventdata,handles)
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
