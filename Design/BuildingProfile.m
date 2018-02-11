function [Cooling, Heating, FanPower,Tzone,Twall,Damper] = BuildingProfile(varargin)
% This function estimates the energy profile of a building (electric, heating and cooling kW)
% Build is a structure of building parameters
% Weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% Date is a vector of points in datenum format at which you are requesting the electric cooling & heating power
Build = varargin{1};
Date = varargin{2};
InternalGains = varargin{3};
ExternalGains = varargin{4};
Tdb = varargin{5};
RH = varargin{6};

% Handle Date input
[s1,~] = size(Date);
if s1 == 1
    Date = Date';
end
nS = length(Date); % number timesteps
Tzone = zeros(nS+1,1);
Twall = zeros(nS+1,1);
Tsupply = zeros(nS,1);
Tzone(1) = varargin{7};
Twall(1) = varargin{8};
TsetC = loadSched(Build,Date,'TsetC');
TsetH = loadSched(Build,Date,'TsetH');
if nS == 1
    dt = 1;
else
    dt1 = Date(2) - Date(1);
    dt = (24*3600) * (Date - [Date(1)-dt1;Date(1:end-1)]); % duration of each time segment [seconds]
end

Cooling = zeros(nS,1);
Heating = zeros(nS,1);
AirFlow = zeros(nS,1);
Damper = zeros(nS,1);

%% Heating/cooling dynamics
%need to step through time to account for moments with no heating/cooling
%as the building moves between the heat setpoint and the cooling setpoint
%note the conversion factor of 5.68 hft^2*F/Btu per 1 m^2*K/W
for t = 1:1:nS  
    if ~Build.VariableStruct.swOffPeakHvac && (TsetC(t) > 26 || TsetH(t) < 19) && (TsetC(max(1,t-round(3600/dt(t)))) > 26 || TsetH(max(1,t-round(3600/dt(t)))) < 19)
        %%HVAC is 'off'
        Tsupply(t) = nan;
        if Tzone(t)<21.5
            Heating(t) = 0.25*3.7382; %constant heating with no air flow (why does Eplus do this?)
        end
    else
        [Heating(t),Cooling(t),Damper(t),AirFlow(t),Tsupply(t)] = HVAClogic(Build,Tzone(t),Twall(t),Tdb(t),RH(t),TsetC(t),TsetH(t),InternalGains(t),ExternalGains(t),dt(t));
    end
    [Heating(t),Cooling(t),Tzone(t+1),Twall(t+1)] = BuildingDynamicResponse(Build,Tdb(t),RH(t),dt(t),InternalGains(t),ExternalGains(t),Cooling(t),Heating(t),Tsupply(t),AirFlow(t),Damper(t),Tzone(t),Twall(t));
end
Cooling(abs(Cooling) < 1e-2) = 0;
FanPower = AirFlow*Build.VariableStruct.FanPower;
end%Ends function BuildingProfile

function [Heating,Cooling,Damper,AirFlow,Tsupply] = HVAClogic(Build,Tzone,Twall,Tdb,RH,TsetC,TsetH,InternalGains,ExternalGains,dt)
rho_Air = 1.225; % kg/m^3
% Fresh = true; %require fresh air recirculation
Damper = 1;
Heating = 0;
Cooling = 0;
Tsupply = nan;
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
MinFlow = Build.VariableStruct.Volume*rho_Air*Build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
Flow = MinFlow;
NetGain = HVAC_equilibrium(Build,Twall,Tdb,TsetH,InternalGains,ExternalGains,dt);
if Build.VariableStruct.swFixFan %small office
    Flow = 10*MinFlow;
    Damper = .1;
end
if (min(TsetH,Build.VariableStruct.ColdAirSet) - Tdb)*Cp_amb*MinFlow*dt - NetGain > 0 
    %% Actively heating  
    if Build.VariableStruct.swFixFan %small office
        Tsupply =  TsetH - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this cooling
    else
        Tsupply = TsetH - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this heating
        Tsupply = min(65,Tsupply); % temperature of supply air to provide this heating 
        Flow = max(MinFlow,NetGain/((TsetH - Tsupply)*Cp_build*dt)); % mass flow of air to provide this heating
        Damper = MinFlow/Flow;
    end
else
    NetGain = HVAC_equilibrium(Build,Twall,Tdb,TsetC,InternalGains,ExternalGains,dt);
    if (Tdb - TsetC)*Cp_amb*MinFlow*dt + NetGain > 0 %cooling if flowing minimum ambient air is insuficient to meet energy removal %net energy beyond what is needed to reach TsetC kJ %%
        %% Actively cooling
        if Build.VariableStruct.swFixFan %small office
            Tsupply =  TsetC - NetGain/(Flow*Cp_amb*dt); % temperature of supply air to provide this cooling
        else
            Tsupply = Build.VariableStruct.ColdAirSet;
            Flow = max(MinFlow, NetGain/((TsetC-Tsupply)*Cp_build*dt)); % mass flow of air to provide this cooling
            if Tdb < Tzone % find economizer position
                Damper = max(MinFlow/Flow,min(1,(Tsupply - TsetC)./(Tdb - TsetC)));
            else
                Damper = MinFlow/Flow; % treat as little ouside air as possible
            end
        end
    else
        %% Passive ventilation with exterior air
        if ~Build.VariableStruct.swFixFan 
            Damper = min(1,(Build.VariableStruct.ColdAirSet - Tzone)/(Tdb - Tzone));%if outdoor air is too cold, dilute with recirculation air
            Flow = MinFlow/Damper;
        end
    end
end
%% De-Humidification
if Build.VariableStruct.swDehumid && m_v_air>m_v_set %estimate de-humidification with Cp (need to add latent heat)
%             Cooling  = Cooling + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
%             Heating  = Heating + Flow*Cp_Air*(Tmix - Build.VariableStruct.DPset)/n; %dehumidification energy in kJ/s
    Latent = Latent_H2O*Flow*Damper*max(0,m_v_air - m_v_set); %latent heat of condensation in kJ/s
    Cooling  = Cooling + Latent;
end
AirFlow = Flow/rho_Air; % flow rate in m^3/s
end%Ends function HVAClogic

function [Heating,Cooling,Tzone,Twall] = BuildingDynamicResponse(Build,Tdb,RH,dt,InternalGains,ExternalGains,Cooling,Heating,Tsupply,AirFlow,Damper,T1,T2)
%%Note this function is very similar to, but cannot be replaced by BuildingSimulate
%%This finds the heating to achieve a Tsupply, where BuildingSimulate
%%already knows the Heating and Cooling
rho_Air = 1.225; % kg/m^3
n = floor(dt/30);% 30 second steps size for smotth simulation of temperature profile
f = 0.1;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = RH/100.*satP; % kPa
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K
%% Specific heat of air
Tdp_set = 273+Build.VariableStruct.DPset;
P_H2O_dp = exp((-5.8002206e3)./Tdp_set + 1.3914993 - 4.8640239e-2*Tdp_set + 4.1764768e-5*Tdp_set.^2 - 1.4452093e-8*Tdp_set.^3 + 6.5459673*log(Tdp_set))/1000; %saturated water vapor pressure at dewpoint
m_v_set = .621945*(P_H2O_dp/(P-P_H2O_dp));%mass fraction of water in air at desired RH
Cp_build = 1.006 + 1.86*m_v_set; %kJ/kg
%%
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
UA_wall = (Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
Tair = zeros(n,1);
Tsolid = zeros(n,1);
Tair(1) = T1;
Tsolid(1) = T2;
AirCapacitance = 2*Build.Area;
Flow = AirFlow*rho_Air; % flow rate in m^3/s
if isnan(Tsupply)
    mode='passive';
else
    mode = 'active';
end
for i = 1:1:n
    Tmix = Damper*Tdb + (1-Damper)*Tair(i);
    Cp_Air = Damper*Cp_amb + (1-Damper)*Cp_build;
    if strcmp(mode,'passive')
        Tsupply = Tmix;
    else
        Q = Cp_Air*(Tsupply-Tmix)*Flow/n;
        if Q>0
            Heating = Heating + Q;
        else
            Cooling = Cooling - Q;
        end
    end
    if Flow ==0
        Tair(i+1) = Tair(i) + ((Tdb - Tair(i))*UA_window + (Tsolid(i) - Tair(i))*(1/f)*UA_wall + InternalGains + Heating)*dt/(n*AirCapacitance); %net change in temperature
    else
        Tair(i+1) = Tair(i) + ((Tdb - Tair(i))*UA_window + (Tsolid(i) - Tair(i))*(1/f)*UA_wall + InternalGains + (Tsupply - Tair(i))*Cp_Air*Flow)*dt/(n*AirCapacitance); %net change in temperature
    end
    Tsolid(i+1) = Tsolid(i) + ((Tdb - Tsolid(i))*(1/(1-f))*UA_wall + (Tair(i+1) - Tsolid(i))*(1/f)*UA_wall + ExternalGains + (Tsky^4 - (273+Tsolid(i))^4)*sigA).*dt/(n*Build.VariableStruct.Capacitance*Build.Area);
end
Tzone = Tair(end);
Twall = Tsolid(end);
end%Ends function BuildingDynamicResponse