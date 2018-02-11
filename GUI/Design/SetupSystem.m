function SetupSystem(str)
% This function updates the GUI graphic to represent the current plant
% configuration, showing buttons and lines connecting generating sources to
% the apropriate bus (AC, DC, thermal...)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
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
    case 'Chiller'
        isType = strcmp('Chiller',type);
    case 'Water Heater'
        isType = strcmp('Boiler',type);
%     Following buttons don't work until we decided what they display
%     case 'Heating Demands'
%         GENINDEX = -1;
%     case 'Hot Water Demands'
%         GENINDEX = -2;
%     case 'Cooling Demands'
%         GENINDEX = -3;
%     case 'AC / DC Conversion'
%         GENINDEX = -4;
end
isType = nonzeros(linspace(1,nG,nG)'.*isType);

% Decide whether to show the user a list of the chosen type
if ~isempty(isType)
    if length(isType) ==1
        GENINDEX = isType;
        EditSystem(handles)
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