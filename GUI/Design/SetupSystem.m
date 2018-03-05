function SetupSystem(str)
% This function updates the GUI graphic to represent the current plant
% configuration, showing buttons and lines connecting generating sources to
% the apropriate bus (AC, DC, thermal...)
global testSystems SYSINDEX GENINDEX
handles = guihandles;
% Put Plant.Generator information into lists
nG = length(testSystems(SYSINDEX).Generator);
isType = zeros(nG,1);
type = cell(nG,1);
name = cell(nG,1);
source = cell(nG,1);
for i = 1:nG
    type(i) = {testSystems(SYSINDEX).Generator(i).Type};
    name(i) = {testSystems(SYSINDEX).Generator(i).Name};
    source(i) = {testSystems(SYSINDEX).Generator(i).Source};
end

% Find indices of relevant components
switch str
    case 'DCgen' 
        isType = strcmp('CHP Generator',type) + strcmp('Electric Generator',type);
        for i = 1:nG
            if isType(i) && ~isfield(testSystems(SYSINDEX).Generator(i).Output,'DirectCurrent')
                isType(i) = false;
            end
        end
    case 'ACgen'
        isType = strcmp('CHP Generator',type) + strcmp('Electric Generator',type);
        for i = 1:nG
            if isType(i) && ~isfield(testSystems(SYSINDEX).Generator(i).Output,'Electricity')
                isType(i) = false;
            end
        end
    case 'Utility'
        isType = strcmp('Utility',type);
    case 'Solar PV'
        isType = strcmp('Solar',type);
    case 'Heater'
        isType = strcmp('Heater',type);
    case 'TES_Hot' 
        isType = strcmp('Thermal Storage',type).*strcmp('Heat',source);
    case 'TES_Cold' 
        isType = strcmp('Thermal Storage',type).*strcmp('Cooling',source);
    case 'Battery'
        isType = strcmp('Electric Storage',type);
    case 'Chiller'
        isType = strcmp('Chiller',type).*strcmp('Elecicity',source);
    case 'Ab Chiller'
        isType = strcmp('Chiller',type).*strcmp('Heat',source);

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