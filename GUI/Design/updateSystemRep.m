function updateSystemRep(hObject, eventdata, handles)
%% System representation (pushbuttons need to appear/disappear if they exist)
global testSystems SYSINDEX
list = {'pushbuttonHeaterInSys';'textAirHeater_H1';'pushbuttonTES1';'textTES1_H1';'pushbuttonHeatingDemands';'pushbuttonICE_mGT';'textICE_mGT_E';'textICE_mGT_H1';'textICE_mGT_H2';'textICE_mGT_H3';'pushbuttonSolarThermalInSys';...
        'textSolarThermal_H2';'pushbuttonWaterHeaterInSys';'textWaterHeater_H2';'pushbuttonTES2';'textTES2_H2';'pushbuttonFuelCell';'textFuelCell_E';'textFuelCell_H1';'textFuelCell_H2';'pushbuttonHotWaterDemands';...
        'pushbuttonAbChillerInSys';'textAbChill_E';'textAbChill_C1';'textAbChill_C2';'pushbuttonChillerInSys';'textChill_E';'textChill_C2';'pushbuttonTES3';'textTES3_C1';'pushbuttonBatteryInSys';...
        'textBattery_E';'pushbuttonCoolingDemands';'pushbuttonWindInSys';'textWind_E';'pushbuttonSolarSterlingInSys';'pushbuttonSolarPVInSys';'textSolarPV_E';'pushbuttonGrid';'textGrid';};
% Hide everything except grid and AC/DC conversion
for i = 1:1:length(list)
    set(handles.(list{i}),'Visible','off')
end
% Show AC/DC conversion by default
set(handles.pushbuttonACDC,'Visible','on')
set(handles.textACDC_E,'Visible','on')

% Selectively show existing component categories
% Does not handle:
%   Absorption Chiller
%   Solar Stirling
%   Solar Thermal
%   TES 1
%   Wind
nG = length(testSystems(SYSINDEX).Generator);
for i = 1:1:nG
    sys = testSystems(SYSINDEX).Generator(i);
    switch sys.Type
        case 'Utility'
            set(handles.pushbuttonGrid,'Visible','on')
            set(handles.textGrid,'Visible','on')
        case 'Heater'
            set(handles.pushbuttonHeaterInSys,'Visible','on')
            set(handles.textAirHeater_H1,'Visible','on')
        case 'Thermal Storage'
            if strcmp(sys.Source,'Heat')
                set(handles.pushbuttonTES2,'Visible','on')
                set(handles.textTES2_H2,'Visible','on')
            elseif strcmp(sys.Source,'Cooling')
                set(handles.pushbuttonTES3,'Visible','on')
                set(handles.textTES3_C1,'Visible','on')
            end
        case 'CHP Generator'
            if sys.VariableStruct.isFuelCell
                set(handles.pushbuttonFuelCell,'Visible','on')
                set(handles.textFuelCell_E,'Visible','on')
                set(handles.textFuelCell_H1,'Visible','on')
                set(handles.textFuelCell_H2,'Visible','on')
            else
                set(handles.pushbuttonICE_mGT,'Visible','on')
                set(handles.textICE_mGT_E,'Visible','on')
                set(handles.textICE_mGT_H1,'Visible','on')
                set(handles.textICE_mGT_H2,'Visible','on')
            end
        case 'Solar'
            set(handles.pushbuttonSolarPVInSys,'Visible','on')
            set(handles.textSolarPV_E,'Visible','on')
        case 'Boiler'
            set(handles.pushbuttonWaterHeaterInSys,'Visible','on')
            set(handles.textWaterHeater_H2,'Visible','on')
        case 'Chiller'
            set(handles.pushbuttonChillerInSys,'Visible','on')
            set(handles.textChill_E,'Visible','on')
            set(handles.textChill_C2,'Visible','on')
        case 'Electric Storage'
            set(handles.pushbuttonBatteryInSys,'Visible','on')
            set(handles.textBattery_E,'Visible','on')
        case 'Electric Generator'
            set(handles.pushbuttonICE_mGT,'Visible','on')
            set(handles.textICE_mGT_E,'Visible','on')
            set(handles.textICE_mGT_H2,'Visible','on')
    end
end

% Check whether hot/cold demands should be shown
if strcmp(get(handles.textICE_mGT_H2,'Visible'),'on') || ...
        strcmp(get(handles.textSolarThermal_H2,'Visible'),'on') || ...
        strcmp(get(handles.textWaterHeater_H2,'Visible'),'on') || ...
        strcmp(get(handles.textTES2_H2,'Visible'),'on') || ...
        strcmp(get(handles.textFuelCell_H2,'Visible'),'on')
    set(handles.pushbuttonHotWaterDemands,'Visible','on')
end
if strcmp(get(handles.textAirHeater_H1,'Visible'),'on') || ...
        strcmp(get(handles.textTES1_H1,'Visible'),'on') || ...
        strcmp(get(handles.textICE_mGT_H1,'Visible'),'on') || ...
        strcmp(get(handles.textFuelCell_H1,'Visible'),'on')
    set(handles.pushbuttonHeatingDemands,'Visible','on')
end
if strcmp(get(handles.textAbChill_C1,'Visible'),'on') || ...
        strcmp(get(handles.textChill_C2,'Visible'),'on') || ...
        strcmp(get(handles.textTES3_C1,'Visible'),'on')
    set(handles.pushbuttonCoolingDemands,'Visible','on')
end