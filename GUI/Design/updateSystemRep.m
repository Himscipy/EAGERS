function updateSystemRep(handles)
%% System representation (pushbuttons need to appear/disappear if they exist)
global testSystems SYSINDEX TestData
list = {'textGrid';'pushbuttonGrid';'textACBus';'textDCBus';'textACDC_E';...
        'pushbuttonACgen';'textACgen_E';'textACgen_H';'pushbuttonHeaterInSys';...
        'textAirHeater_H';'pushbuttonTES_Hot';'textTES_Hot_H';'pushbuttonHeatingDemands';...
        'pushbuttonDCgen';'textDCgen_E';'textDCgen_H';...
        'pushbuttonAbChill_InSys';'textABchill_H';'textAbChill_C';...
        'pushbuttonChillerInSys';'textChill_E';'textChill_C';...
        'pushbuttonTES_Cold';'textTES_Cold_C';'pushbuttonCoolingDemands';...
        'pushbuttonBatteryInSys';'textBattery_E';...
        'pushbuttonWindInSys';'textWind_E';...
        'textSolarStering_E';'pushbuttonSolarSterlingInSys';...
        'pushbuttonSolarPVInSys';'textSolarPV_E';};
% Hide everything except grid and AC/DC conversion
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','off')
end

% Selectively show existing component categories
% Does not handle:
%   Absorption Chiller
%   Solar Stirling
%   Solar Thermal
%   TES 1
%   Wind
if ~isfield(testSystems(SYSINDEX),'Building')
    testSystems(SYSINDEX).Building = [];
end
if isfield(TestData,'Demand')
    demand_types = fieldnames(TestData.Demand);
else 
    demand_types = {};
end
[testSystems(SYSINDEX).Generator,nG] = check_ac_dc(testSystems(SYSINDEX).Generator,testSystems(SYSINDEX).Building,demand_types);
for i = 1:1:nG
    sys = testSystems(SYSINDEX).Generator(i);
    switch sys.Type
        case 'Utility'
            set(handles.pushbuttonGrid,'Visible','on')
            set(handles.textGrid,'Visible','on')
            set(handles.textACBus,'Visible','on')
        case 'AC_DC'
            set(handles.pushbuttonACDC,'Visible','on')
            set(handles.textACDC_E,'Visible','on')
            set(handles.textACBus,'Visible','on')
            set(handles.textDCBus,'Visible','on')
        case 'Heater'
            set(handles.pushbuttonHeaterInSys,'Visible','on')
            set(handles.textAirHeater_H,'Visible','on')
            set(handles.pushbuttonHeatingDemands,'Visible','on')
        case 'Thermal Storage'
            if strcmp(sys.Source,'Heat')
                set(handles.pushbuttonTES_Hot,'Visible','on')
                set(handles.textTES_Hot_H,'Visible','on')
                set(handles.pushbuttonHeatingDemands,'Visible','on')
            elseif strcmp(sys.Source,'Cooling')
                set(handles.pushbuttonTES_Cold,'Visible','on')
                set(handles.textTES_Cold_C,'Visible','on')
                set(handles.pushbuttonCoolingDemands,'Visible','on')
            end
        case {'CHP Generator';'Electric Generator'}
            if isfield(sys.Output,'DirectCurrent')
                set(handles.pushbuttonDCgen,'Visible','on')
                set(handles.textDCgen_E,'Visible','on')
                set(handles.textDCBus,'Visible','on')
                if isfield(sys.Output,'Heat')
                    set(handles.textDCgen_H,'Visible','on')
                    set(handles.pushbuttonHeatingDemands,'Visible','on')
                end
            else
                set(handles.pushbuttonACgen,'Visible','on')
                set(handles.textACgen_E,'Visible','on')
                set(handles.textACBus,'Visible','on')
                if isfield(sys.Output,'Heat')
                    set(handles.textACgen_H,'Visible','on')
                    set(handles.pushbuttonHeatingDemands,'Visible','on')
                end
            end
        case 'Solar'
            %need to make it on AC or DC bus
            set(handles.pushbuttonSolarPVInSys,'Visible','on')
            set(handles.textSolarPV_E,'Visible','on')
        case 'Chiller'
            if strcmp(sys.Source,'Heat')
                set(handles.pushbuttonAbChill_InSys,'Visible','on')
                set(handles.textABchill_H,'Visible','on')
                set(handles.textABchill_C,'Visible','on')
            elseif strcmp(sys.Source,'Electricity')
                set(handles.pushbuttonChillerInSys,'Visible','on')
                set(handles.textChill_E,'Visible','on')
                set(handles.textChill_C,'Visible','on')
            end
            set(handles.pushbuttonCoolingDemands,'Visible','on')
        case 'Electric Storage'
            %need to make it on AC or DC bus
            set(handles.pushbuttonBatteryInSys,'Visible','on')
            set(handles.textBattery_E,'Visible','on')
    end
end
end%ends function updateSystemRep