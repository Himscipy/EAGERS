function [Equipment, InteriorLighting, ExteriorLighting, Cooling, Heating, FanPower, OtherLoads,Tzone,Twall] = BuildingProfile(varargin)
% This function estimates the energy profile of a building (electric, heating and cooling kW)
% Build is a structure of building parameters
% Weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% Date is a vector of points in datenum format at which you are requesting the electric cooling & heating power

%% Constants
rho_Air = 1.225; % kg/m^3
ExtLights_B = .04;%.125; % fraction of day exterior lights remain on after sunrise

%% Initialize
% Inputs
Build = varargin{1};
Weather = varargin{2};
Date = varargin{3};
if length(varargin) == 4 && isempty(varargin{4})
    Location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
else Location = varargin{4};
end

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


% Schedules
sched = fieldnames(Build.Schedule);
for i = 1:1:length(sched)
    Profile.(sched{i}) = zeros(nS,1);
end
Tdb = zeros(nS,1); % Drybulb temperature
RH = zeros(nS,1); % Humidity ratio
Cooling = zeros(nS,1);
Heating = zeros(nS,1);
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
sgWalls = Build.VariableStruct.WallAbsorption*(1-Build.VariableStruct.WindowWallRatio)*(sgDirect+sgDiffuse)/1000;%solar gain (heat) absorbed by walls
vlWindows = Build.VariableStruct.LightTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%visible light transmitted through windows
sgRoof = Build.VariableStruct.RoofArea*Build.VariableStruct.WallAbsorption*irrad.*max(0,dot(n*[0,0,1],DN,2))/1000;%solar gain (heat) through roof

%% Interior lighting load values
InteriorLighting = Build.Area*Build.VariableStruct.InteriorLights/1000*Profile.interiorlights.*(1 - Build.VariableStruct.DaylightPercent*vlWindows/max(vlWindows)); % kW of lighting load

%% Internal gains
internalGainOccupants = Build.Area * Build.VariableStruct.occupancy * Profile.occupancy * 0.120; % heat from occupants (120 W)
InternalGains = internalGainOccupants + Equipment + InteriorLighting + sgWindows;

%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = RH/100.*satP; % kPa
Tdp = 6.54 + 14.526*log(P_H2O) + 0.7389*log(P_H2O).^2 + 0.09486*log(P_H2O).^3 + 0.4569*(P_H2O).^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K
m_v_air = .621945*(P_H2O./(P-P_H2O));%mass fraction of water in air 

%% Specific heat of air
Tdp_set = 273+Build.VariableStruct.DPset;
P_H2O_dp = exp((-5.8002206e3)./Tdp_set + 1.3914993 - 4.8640239e-2*Tdp_set + 4.1764768e-5*Tdp_set.^2 - 1.4452093e-8*Tdp_set.^3 + 6.5459673*log(Tdp_set))/1000; %saturated water vapor pressure at dewpoint
m_v_set = .621945*(P_H2O_dp/(P-P_H2O_dp));%mass fraction of water in air at desired RH
Cp_build = 1.006 + 1.86*m_v_set; %kJ/kg
Latent_H2O = 2500.8 - 2.36*Build.VariableStruct.DPset +.0016*Build.VariableStruct.DPset^2 - 0.00006*Build.VariableStruct.DPset^3;% kJ/kg of latent heat

% Minimum air flow rate
MinFlow = Build.VariableStruct.Volume*rho_Air*Build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
%% Heating/cooling dynamics
%need to step through time to account for moments with no heating/cooling
%as the building moves between the heat setpoint and the cooling setpoint
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
UA_wall = (Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
%note the conversion factor of 5.68 hft^2*F/Btu per 1 m^2*K/W
Tzone = zeros(nS+1,1);
Twall = zeros(nS+1,1);
Tsupply = zeros(nS,1);

%% assuming starting jan1 need to modify for non full-year simulation
Tzone(1) = Profile.TsetH(1) + 1.4;
Twall(1) = Tzone(1);
mode = 'heating';
%%%%
n = 10;
f = .1;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
Tair = ones(n+1,1);
Tsolid = ones(n+1,1);
AirCapacitance = 2*Build.Area;
for t = 1:1:nS    
    Fresh = true; %require fresh air recirculation
    Damper = 1;
    Flow = MinFlow;
    
    %%corrections for when HVAC is off or not bringing in fresh air
    if ~Build.VariableStruct.swOffPeakHvac
        if (Profile.TsetC(t) > 26 || Profile.TsetH(t) < 19) && (Profile.TsetC(max(1,t-round(3600/dt(t)))) > 26 || Profile.TsetH(max(1,t-round(3600/dt(t)))) < 19)
            Fresh = false;
            Flow = 0; 
            Tsupply(t) = Tzone(t);
            Tzone(t+1) = Tzone(t);
            mode = 'passive';
        end 
    end    
    if Fresh %%need to figure out what is happening when heating or cooling and not fresh.
        if Tzone(t)>Twall(t)
            Tzone(t+1) = max(Twall(t),Profile.TsetH(t));
        else
            Tzone(t+1) = Profile.TsetH(t);
        end
        Twall(t+1) = Twall(t) + ((Tdb(t) - Twall(t))*(1/(1-f))*UA_wall + (Tzone(t+1) - Twall(t))*(1/f)*UA_wall + sgWalls(t) + sgRoof(t) + (Tsky^4 - (273+Twall(t))^4)*sigA).*dt(t)/(Build.VariableStruct.Capacitance*Build.Area);%estimate of wall temperature at t+1
        Twall(t+1) = 0.5*Twall(t) + 0.25*Twall(t+1) + 0.25*(Twall(t) + ((Tdb(t) - Twall(t+1))*(1/(1-f))*UA_wall + (Tzone(t+1) - Twall(t+1))*(1/f)*UA_wall + sgWalls(t) + sgRoof(t) + (Tsky^4 - (273+Twall(t+1))^4)*sigA).*dt(t)/(Build.VariableStruct.Capacitance*Build.Area));%average wall temp now with wall temp at t+1
        NetGain = ((Tdb(t) - Tzone(t+1))*UA_window + (Twall(t+1) - Tzone(t+1))*(1/f)*UA_wall + InternalGains(t))*dt(t); %net energy gain into the zone in kJ, assuming the zone goes to TsetH
        if (min(Profile.TsetH(t),Build.VariableStruct.ColdAirSet) - Tdb(t))*Cp_amb(t)*Flow*dt(t) - NetGain > 0 %% Actively heating  
            mode = 'heating';
            Tsupply(t) = Tzone(t+1) - NetGain/(Flow*Cp_amb(t)*dt(t)); % temperature of supply air to provide this heating
            Tsupply(t) = min(65,Tsupply(t)); % temperature of supply air to provide this heating 
            Flow = max(MinFlow,NetGain/((Tzone(t+1) - Tsupply(t))*Cp_build*dt(t))); % mass flow of air to provide this heating
            Damper = MinFlow/Flow;
        else
            Tzone(t+1) = Profile.TsetC(t);
            NetGain = ((Tdb(t) - Tzone(t+1))*UA_window + (Twall(t+1) - Tzone(t+1))*(1/f)*UA_wall + InternalGains(t))*dt(t); %net energy gain into the zone in kJ, assuming the zone goes to TsetC
            if (Tdb(t) - Profile.TsetC(t))*Cp_amb(t)*Flow*dt(t) + NetGain > 0 %cooling if flowing minimum ambient air is insuficient to meet energy removal %net energy beyond what is needed to reach TsetC kJ %%
                mode = 'cooling';%Actively cooling
                Tsupply(t) = Build.VariableStruct.ColdAirSet;
                Flow = max(MinFlow, NetGain/((Tzone(t+1)-Tsupply(t))*Cp_build*dt(t))); % mass flow of air to provide this cooling
                if Tdb(t) < Tzone(t) % find economizer position
                    Damper = max(MinFlow/Flow,min(1,(Tsupply(t) - Tzone(t+1))./(Tdb(t) - Tzone(t+1))));
                else
                    Damper = MinFlow/Flow; % treat as little ouside air as possible
                end
            else%%Passive ventilation with exterior air
                Tzone(t+1) = Tzone(t);
                NetGain = ((Tdb(t) - Tzone(t+1))*UA_window + (Twall(t+1) - Tzone(t+1))*(1/f)*UA_wall + InternalGains(t))*dt(t); %net energy gain into the zone in kJ, assuming the zone stays the same
                mode = 'passive';
                Damper = min(1,(Build.VariableStruct.ColdAirSet - Tzone(t))/(Tdb(t) - Tzone(t)));%if outdoor air is too cold, dilute with recirculation air
                Flow = MinFlow/Damper;
            end
        end
   
    end
    
    if Build.VariableStruct.swFixFan %small office
        if Fresh
            Flow = 10*MinFlow;
            Damper = .1;
            Tsupply(t) =  Tzone(t+1) - NetGain/(Flow*Cp_amb(t)*dt(t)); % temperature of supply air to provide this cooling
        else
            Flow = 0; 
            Tsupply(t) = Tzone(t);
            if Tzone(t+1)<21.5
                Heating(t) = 0.25*3.7382; %constant heating with no air flow (why does Eplus do this?)
            end
        end 
    end
    DirectHeating = Heating(t);
    Tair(1) = Tzone(t);
    Tsolid(1) = Twall(t);
    for i = 1:1:n
        Tmix = Damper*Tdb(t) + (1-Damper)*Tair(i);
        Cp_Air = Damper*Cp_amb(t) + (1-Damper)*Cp_build;
        if strcmp(mode,'passive')
            Tsupply(t) = Tmix;
        else
            Q = Cp_Air*(Tsupply(t)-Tmix)*Flow/n;
            if Q>0
                Heating(t) = Heating(t) + Q;
            else
                Cooling(t) = Cooling(t) - Q;
            end
        end
        
        %% De-Humidification
        if Build.VariableStruct.swDehumid && m_v_air(t)>m_v_set %estimate de-humidification with Cp (need to add latent heat)
%             Cooling(t)  = Cooling(t) + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
%             Heating(t)  = Heating(t) + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
            Latent = 0.15*Latent_H2O*Flow*Damper*max(0,m_v_air(t) - m_v_set); %latent heat of condensation in kJ/s
            Cooling(t)  = Cooling(t) + Latent;
        end
        
        
        Tair(i+1) = Tair(i) + ((Tdb(t) - Tair(i))*UA_window + (Tsolid(i) - Tair(i))*(1/f)*UA_wall + InternalGains(t) + (Tsupply(t) - Tair(i))*Cp_Air*Flow + DirectHeating)*dt(t)/(n*AirCapacitance); %net change in temperature
        %%correct temperature & heating value
        if t>1
            if Profile.TsetC(t-1)<=Profile.TsetC(t) && Tair(i+1)>Profile.TsetC(t)
                Q = (Tair(i+1) - Profile.TsetC(t))*AirCapacitance/dt(t);
%                 Tsupply(t) = Tsupply(t) - Q/(Cp_Air*Flow);
                Tair(i+1) = Profile.TsetC(t);
                if Heating(t)>0 
                    dQ = min(Heating(t),Q);
                    Heating(t) = Heating(t) - dQ;
                else dQ = 0;
                end
                Cooling(t) = Cooling(t) + (Q-dQ);
            elseif Profile.TsetH(t-1)>= Profile.TsetH(t) && Tair(i+1)<Profile.TsetH(t)
                Q = (Profile.TsetH(t) - Tair(i+1))*AirCapacitance/dt(t);
%                 Tsupply(t) = Tsupply(t) + Q/(Cp_Air*Flow);
                Tair(i+1) = Profile.TsetH(t);
                if Cooling(t)>0 
                    dQ = min(Cooling(t),Q);
                    Cooling(t) = Cooling(t) - dQ;
                else dQ = 0;
                end
                Heating(t) = Heating(t) + (Q-dQ);
            end
        end
        Tsolid(i+1) = Tsolid(i) + ((Tdb(t) - Tsolid(i))*(1/(1-f))*UA_wall + (Tair(i+1) - Tsolid(i))*(1/f)*UA_wall + sgWalls(t) + sgRoof(t) + (Tsky^4 - (273+Tsolid(i))^4)*sigA).*dt(t)/(n*Build.VariableStruct.Capacitance*Build.Area);
    end
    Tzone(t+1) = Tair(end);
    Twall(t+1) = Tsolid(end);
    AirFlow(t) = Flow/rho_Air; % flow rate in m^3/s

%     %%full enthalpy calculations
%     if Tdp(t)>Build.VariableStruct.DPset %must dehumidify incoming air
%         %ambient air
%         AmbientAir = makeAir(Tdb(t),RH(t),Damper*Flow,'rel');
%         %recirculated air
%         RecircAir = makeAir(Tact,Build.VariableStruct.DPset,(1-Damper)*Flow,'dp');
%         %mixed air
%         MixedAir = MixAir(RecircAir,AmbientAir);
%         CooledAir = MixedAir;
%         CooledAir.T = Build.VariableStruct.DPset+273.15;
%         sat_atDP = makeAir(Build.VariableStruct.DPset,100,Flow,'rel');
%         CooledAir.H2O = sat_atDP.H2O; %remove water
%         Cooling(t)  = (-enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) + enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18)); %dehumidification energy in kJ/s
%         HeatedAir = CooledAir;
%         HeatedAir.T = Tsupply(t)+273.15;
%         Heating(t)  = (enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18)); %reheat energy in kJ/s
%     end
    
end
Tzone = Tzone(2:end);
Twall = Twall(2:end);
Cooling(abs(Cooling) < 1e-2) = 0;
FanPower = AirFlow*Build.VariableStruct.FanPower;

if isfield(Build.VariableStruct,'DataCenter')%data center is an equipment load but not an internal gain into the zones
    Equipment = Equipment + Build.VariableStruct.DataCenter*Profile.datacenter;
end
end%Ends function BuildingProfile


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
