function [Tzone,Equipment, InteriorLighting, ExteriorLighting, OtherLoads] = BuildingSimulate(varargin)
% This function estimates the temperature of a building given the ambient conditions, heating and cooling
% Build is a structure of building parameters
% Weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% Date is a vector of points in datenum format at which you are requesting the electric cooling & heating power
% Location is a structure with the latitude, longitude, and time zone
% Cooling is the chiller electric power in kW
% Heating is the heater power in kW
% Tact is the initial zone temperature
%% Constants
rho_Air = 1.225; % kg/m^3
ExtLights_B = .04;%.125; % fraction of day exterior lights remain on after sunrise

%% Initialize
% Inputs
Build = varargin{1};
Weather = varargin{2};
Date = varargin{3};
Location = varargin{4};
Cooling = varargin{5};
Heating = varargin{6};
FanPower = varargin{7};
Tact = varargin{8};
if length(varargin)>8
    Tbuild = varargin{9};
end
Heating = Heating*Build.VariableStruct.COP_H;%convert from kWe to kW thermal
Cooling = Cooling*Build.VariableStruct.COP_C;%convert from kWe to kW thermal

AirFlow = FanPower/Build.VariableStruct.FanPower;
Flow = AirFlow*rho_Air;
% Handle Date input
[s1,~] = size(Date);
if s1 == 1
    Date = Date';
end

% Time
[Y,Month,D,H,M,S] = datevec(Date);
h1 = H(1) + M(1) + S(1);
D(H==0&M==0&S==0) = D(H==0&M==0&S==0) - 1;
H(H==0&M==0&S==0) = 24; % change hour 0 to hour 24 for interpolating.
if h1 == 0 % put start date back.
    D(1) = D(1) + 1;
    H(1) = 0;
end
days = ceil(Date(end) - datenum([Y(1),Month(1),D(1)]));
daysAfter1_1_17 = floor(Date(1)) - datenum([2017,1,1]);
wd = 1 + mod(daysAfter1_1_17,7); % day of the week, sunday is 1, saturday is 7
sd = floor(Date(1)); % start date
nS = length(Date); % number timesteps
dt1 = Date(2) - Date(1);
dt = (24*3600) * (Date - [Date(1)-dt1;Date(1:end-1)]); % duration of each time segment [seconds]

% Minimum air flow rate
MinFlow = Build.VariableStruct.Volume*rho_Air*Build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow


% Schedules
sched = fieldnames(Build.Schedule);
for i = 1:1:length(sched)
    Profile.(sched{i}) = zeros(nS,1);
end
Tdb = zeros(nS,1); % Drybulb temperature
RH = zeros(nS,1); % Humidity ratio
AirFlow = zeros(nS,1);

%% Fill out schedule values for each timestep in input Date
% If next point in schedule is shorter than ramp, then ramp shortens to half the gap
h8761 = (0:1:8760)';
p = 1;
DateIndex = cell(days);
for d = 1:1:days
    h_of_y = 24*(Date(p) - datenum([Y(p),1,1])); % hour since start of year
    % Load schedules for this day
    Out = loadSched(Build,wd,h_of_y);
    p2 = p;
    while p2<nS && Date(p2)<=sd+d
        p2 = p2+1;
    end
    if Date(p2)>sd+d
        p2 = p2-1;
    end
    for i = 1:1:length(sched)
        Profile.(sched{i})(p:p2) = interp1(Out.(sched{i})(:,1), ...
            Out.(sched{i})(:,2), H(p:p2)+M(p:p2)/60+S(p:p2)/3600);
    end
    h_of_y = 24*(Date(p:p2) - datenum([Y(p),1,1])); % hour since start of year
    Tdb(p:p2) = interp1(h8761,[Weather.Tdb(1);Weather.Tdb],h_of_y);
    RH(p:p2) = interp1(h8761,[Weather.RH(1);Weather.RH],h_of_y);
    DateIndex(d) = {p:p2};
    p = p2+1;
    wd = wd+1;
    if wd == 8
        wd = 1;
    end
end

%% Equipment load values
Equipment = Build.Area*Build.VariableStruct.equipment/1000*Profile.equipment; % kW of equipment load

%% Other load values
OtherLoads = zeros(nS,1);
if isfield(Build.VariableStruct,'WaterSystem')
    OtherLoads = Build.Area*Build.VariableStruct.WaterSystem*Profile.watersystem;
end

%% Exterior lighting load values
if isfield(Profile,'exteriorlights_solarcontroled')
    for day = 1:days
        currDate = floor(Date(DateIndex{day}(1)));
        HoD = 24*(Date(DateIndex{day}) - currDate); % hour of day
        frac_dark = ones(length(HoD),1);
        [Sunrise, Sunset, ~, ~] = SolarCalc(Location.Longitude, ...
            Location.Latitude, Location.TimeZone, currDate);
        
        % Find hour at which sunrise occurs
        h_sr = 1;
        while HoD(h_sr)<(Sunrise*24+ExtLights_B) && h_sr<length(HoD)
            h_sr = h_sr+1;
        end
        if h_sr == 1
            frac_dark(1) = (Sunrise*24+ExtLights_B)/HoD(1);
        else
            frac_dark(h_sr) = ((Sunrise*24+ExtLights_B) - HoD(h_sr-1))/(HoD(h_sr)-HoD(h_sr-1));
        end
        
        % Find hour at which sunset occurs
        h_ss = h_sr+1;
        if h_ss < length(HoD)
            while HoD(h_ss)<(Sunset*24-ExtLights_B)  && h_ss<length(HoD)
                h_ss = h_ss+1;
            end
            frac_dark(h_ss) = 1-((Sunset*24-ExtLights_B) - HoD(h_ss-1))/(HoD(h_ss)-HoD(h_ss-1));
        end
        frac_dark(h_sr+1:h_ss-1) = 0;
        Profile.exteriorlights_solarcontroled(DateIndex{day}) = min(Profile.exteriorlights_solarcontroled(DateIndex{day}'),frac_dark);
    end
    ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights_fixed + Build.VariableStruct.ExteriorLights*Profile.exteriorlights_solarcontroled;
else
    ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights;
end


%% Solar gain
[sY, ~, ~, ~, ~, ~] = datevec(Date(1));
sD = datenum([sY,1,1,0,0,0]);
nD = length(Date);
or = Build.VariableStruct.Orientation;
irrad = interp1(linspace(0,8760,8761),[Weather.irradDireNorm(end);Weather.irradDireNorm],mod((Date-sD)*24,8760));
sgDiffuse = Build.VariableStruct.WallArea*interp1(linspace(0,8760,8761),[Weather.irradDiffHorz(end);Weather.irradDiffHorz],mod((Date-sD)*24,8760));%Diffuse irradiance (W)
[~, ~, azimuth, zenith] = SolarCalc(Location.Longitude, Location.Latitude, Location.TimeZone, Date);
DN = [cosd(azimuth).*cosd(90 - zenith), sind(azimuth).*cosd(90 - zenith),sind(90 - zenith)];%Direct normal (incoming vector of sunlight)
n = ones(nD,1);
sgDirect = 0.25*Build.VariableStruct.WallArea*irrad.*(max(0,dot(n*[cosd(or),sind(or),0],DN,2)) + max(0,dot(n*[-cosd(or),-sind(or),0],DN,2)) + max(0,dot(n*[-sind(or),cosd(or),0],DN,2)) + max(0,dot(n*[sind(or),-cosd(or),0],DN,2)));
sgWindows = Build.VariableStruct.WindowTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%solar gain (heat) through windows http://www.commercialwindows.org/vt.php
sgWalls = Build.VariableStruct.WallTransmittance*(1-Build.VariableStruct.WindowWallRatio)*(sgDirect+sgDiffuse)/1000;%solar gain (heat) through walls
vlWindows = Build.VariableStruct.LightTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%visible light transmitted through windows
sgRoof = Build.VariableStruct.RoofArea*Build.VariableStruct.WallTransmittance*irrad.*max(0,dot(n*[0,0,1],DN,2))/1000;%solar gain (heat) through roof

%% Interior lighting load values
InteriorLighting = Build.Area*Build.VariableStruct.InteriorLights/1000*Profile.interiorlights.*(1 - Build.VariableStruct.DaylightPercent*vlWindows/max(vlWindows)); % kW of lighting load

%% Internal gains
internalGainOccupants = Build.Area * Build.VariableStruct.occupancy * Profile.occupancy * 0.120; % heat from occupants (120 W)
InternalGains = internalGainOccupants + Equipment + InteriorLighting + (sgWindows + sgWalls + sgRoof);

%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = RH/100.*satP; % kPa
Tdp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K

%% Specific heat of air
% Tsupply = 22+273.15;
% P = 101.325;
% satP = exp((-5.8002206e3)./Tset + 1.3914993 - 4.8640239e-2*Tset + 4.1764768e-5*Tset.^2 - 1.4452093e-8*Tset.^3 + 6.5459673*log(Tset))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
% P_H2O = exp(5423*(1/273-1./(273+Build.VariableStruct.DPset)))./exp(5423*(1/273-1./(273+Tset))).*satP;%Clausius-Clapeyron equation to calculate relative humidity from dewpoint temperature
% w = (P_H2O/(P-P_H2O));
w = 0.0085;
Cp_build = 1.006 + 1.86*(.621945*w); %kJ/kg

%% Heating/cooling dynamics
%need to step through time to account for moments with no heating/cooling
%as the building moves between the heat setpoint and the cooling setpoint
HeatTransmittance = (Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue + Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
%note the conversion factor of 5.68 hft^2*F/Btu per 1 m^2*K/W
Tzone = zeros(nS,1);
n = 10;
for t = 1:1:nS
    
    
%     %% De-Humidification
%     if Build.VariableStruct.swDehumid
%         %estimate de-humidification with Cp
%         if Tdp(t)>Build.VariableStruct.DPset %must dehumidify incoming air
%             Cp_Air = mean(Damper*Cp_amb(t) + (1-Damper)*Cp_build);
%             Cooling(t)  = max(0,Cooling(t) - Flow(t)*Cp_amb*(Tmix - Build.VariableStruct.DPset)); %dehumidification energy in kJ/s
%        end
%     end
    Tprof = linspace(Tact,Profile.TsetH(t),n)';
    Energy2Add = sum(-(Tdb(t) - Tprof)*HeatTransmittance - InternalGains(t))*dt(t)/n - (Tact - Profile.TsetH(t))*Build.VariableStruct.Capacitance*Build.Area; %net energy addition needed in kJ
    if Energy2Add - sum((Tdb(t) - Tprof)*Cp_amb(t)*MinFlow)*dt(t)/n> 0 %% Actively heating if zone + minimum fresh air needs heat
        Flow(t) = max(MinFlow,Energy2Add/((Tsupply - mean(Tprof))*Cp_build*dt(t))); % mass flow of air to provide this heating
        Damper = MinFlow/Flow(t);
        Tmix = mean(Damper*Tdb(t) + (1-Damper).*Tact);
        Cp_Air = mean(Damper*Cp_amb(t) + (1-Damper)*Cp_build);
        Tsupply = Heating(t)/(Cp_Air*Flow(t))+Tmix;
    else
        Tprof = linspace(Tact,Profile.TsetC(t),n)';
        Energy2Remove = sum((Tdb(t) - Tprof)*HeatTransmittance + InternalGains(t))*dt(t)/n + (Tact - Profile.TsetC(t))*Build.VariableStruct.Capacitance*Build.Area; %net energy needed to be removed in kJ
        if Energy2Remove>0
            Tsupply = Build.VariableStruct.ColdAirSet;
            Flow(t) = max(MinFlow,Energy2Remove/((mean(Tprof)-Tsupply)*Cp_build*dt(t))); % mass flow of air to provide this cooling
            Tmix = Tsupply + Cooling(t)/(Cp_Air*Flow(t));
            if Tdb(t) < Tact % find economizer position
                Damper = max(min(1,MinFlow/Flow(t)),min(1,(Build.VariableStruct.ColdAirSet - Tact)./(Tdb(t) - Tact)));
            else
                Damper = MinFlow/Flow(t); % treat as little ouside air as possible
            end
        else %passive
            Damper = min(1,(Build.VariableStruct.ColdAirSet - Tact)/(Tdb(t) - Tact));%if outdoor air is too cold, dilute with recirculation air
            Cp_Air = mean(Damper*Cp_amb(t) + (1-Damper)*Cp_build);
%             Flow(t) = MinFlow/Damper;
            Tmix = mean(Damper*Tdb(t) + (1-Damper).*Tact);
            if Flow(t) == 0
                Tsupply = Tmix;
            else
                Tsupply = (Heating(t)+Cooling(t))/(Cp_Air*Flow(t))+Tmix;
            end
        end
    end
    for i = 1:1:n
        Tact = Tact + ((Tdb(t) - Tact)*HeatTransmittance + InternalGains(t) + (Tsupply - Tact)*Cp_Air*Flow(t)).*dt(t)/(n*Build.VariableStruct.Capacitance*Build.Area); %net change in temperature
    end
    if ~Build.VariableStruct.swOffPeakHvac
        if Profile.TsetC(t) > 26 || Profile.TsetH(t) < 19
           if Tact < Profile.TsetH(t)
               Tact = Profile.TsetH(t);
           end
           if Tact > Profile.TsetC(t)
               Tact = Profile.TsetC(t);
           end
        end
    end
    Tzone(t) = Tact;
end
if length(varargin)>8
    figure
    plot(Date-Date(1),Tbuild,'b')
    hold on
    plot(Date-Date(1),Tzone,'r')
end

if isfield(Build.VariableStruct,'DataCenter')%don't include data center in internal gains
    Equipment = Equipment + Build.VariableStruct.DataCenter*Profile.datacenter;
end

end%Ends function BuildingSimulate


function Out = loadSched(Build,wd,h_of_y)
sched = fieldnames(Build.Schedule);
for i = 1:1:length(sched)
    s = nnz(Build.Schedule.(sched{i}).Seasons<=h_of_y)+1;% get season
    if isfield(Build.Schedule.(sched{i}),'Ramp')
        Ramp = Build.Schedule.(sched{i}).Ramp;
    else
        Ramp = 1e-4;
    end
    % get schedule
    if wd == 1% Sunday
        day = 'Sun';
    elseif wd == 7 % Saturday
        day = 'Sat';
    else % Weekday
        day = 'Weekday';
    end
    if iscell(Build.Schedule.(sched{i}).(day))
        Out.(sched{i}) = convertSched(Build.Schedule.(sched{i}).(day){s},Ramp);
    else
        Out.(sched{i}) = convertSched(Build.Schedule.(sched{i}).(day),Ramp);
    end
end

end%Ends function loadSched


function newsched = convertSched(sched,Ramp)
[m,n] = size(sched);
if m==2
    newsched = sched; %this is constant all day, already made 0 hour and 24 hour
else
    newsched = zeros(2*(m-1),n);
    sched(1,2) = sched(2,2);
    newsched(1,:) = sched(1,:);%hour 0
    newsched(end,:) = sched(m,:);%hour 24
    if Ramp<1e-3
        for i = 2:1:m-1 
            newsched(2*i-2,1) = sched(i);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+Ramp;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    else 
        for i = 2:1:m-1 %add points in the middle so it can be properly interpolated
            t_bef = sched(i) - sched(i-1);
            t_aft = sched(i+1)-sched(i);
            Ramp2 = min([Ramp, t_aft/2,t_bef/2]);
            newsched(2*i-2,1) = sched(i)-0.5*Ramp2;
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+0.5*Ramp2;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    end
end

end%Ends function convertSched
