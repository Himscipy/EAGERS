function [Tzone,Twall] = building_simulate(Build,Tdb,RH,dt,InternalGains,ExternalGains,Cooling,Heating,AirFlow,Damper,Tzone_init,Twall_init)
% This function estimates the temperature of a building given the ambient conditions, heating and cooling
% Build is a structure of building parameters
% Cooling is the chiller electric power in kW
% Heating is the heater power in kW
nS = length(dt); % number timesteps
Tzone = zeros(nS+1,1);
Twall = zeros(nS+1,1);
Tzone(1) = Tzone_init;
Twall(1) = Twall_init;

%% Ambient dewpoint
rho_Air = 1.225; % kg/m^3
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

%% Heating/cooling dynamics
%need to step through time to account for moments with no heating/cooling
%as the building moves between the heat setpoint and the cooling setpoint
UA_window = Build.VariableStruct.WallArea*Build.VariableStruct.WindowWallRatio*Build.VariableStruct.WindowUvalue/1000;
UA_wall = (Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio)/Build.VariableStruct.WallRvalue + Build.VariableStruct.RoofArea/Build.VariableStruct.RoofRvalue)/1000;% Transmittance in kW/K, see the following about whole wall R-Value http://web.ornl.gov/sci/buildings/docs/Thermal-Performance-and-Wall-Ratings.pdf
sigA = Build.VariableStruct.WallEmissivity*(Build.VariableStruct.WallArea*(1-Build.VariableStruct.WindowWallRatio) + Build.VariableStruct.RoofArea)*5.670367e-8/1000;% kW/K^4
Tsky = (16+273); %temperature in K that building is exhanging radiative heat transfer with
%note the conversion factor of 5.68 ft^2*F/Btu per 1 m^2*K/W

f = .1;%fraction of wall R-value attributed to between wall and zone, remainder of R-value is between ambient and wall (needs to be 0.1 for small office)
AirCapacitance = 2*Build.Area;
for t = 1:1:nS
    Flow = AirFlow(t)*rho_Air; % flow rate in m^3/s
    if Build.VariableStruct.swDehumid && m_v_air(t)>m_v_set %estimate de-humidification with Cp (need to add latent heat)
        Latent = Latent_H2O*Flow*Damper(t)*max(0,m_v_air(t) - m_v_set); %latent heat of condensation in kJ/s
        Cooling(t)  = Cooling(t) - Latent;
    end
    n = floor(dt(t)/30);% 30 second steps size for smotth simulation of temperature profile
    Tair = ones(n+1,1);
    Tsolid = ones(n+1,1);
    Tair(1) = Tzone(t);
    Tsolid(1) = Twall(t);
    for i = 1:1:n
        if Flow == 0
            %Direct heating mode (EnergyPlus adds heat without any HVAC flow)
            Tair(i+1) = Tair(i) + ((Tdb(t) - Tair(i))*UA_window + (Tsolid(i) - Tair(i))*(1/f)*UA_wall + InternalGains(t) + Heating(t))*dt(t)/(n*AirCapacitance); %net change in temperature
        else
            Tmix = Damper(t)*Tdb(t) + (1-Damper(t))*Tair(i);
            Cp_Air = Damper(t)*Cp_amb(t) + (1-Damper(t))*Cp_build;
            Tsupply = (Heating(t)-Cooling(t))/(Flow*Cp_Air) + Tmix;
            Tair(i+1) = Tair(i) + ((Tdb(t) - Tair(i))*UA_window + (Tsolid(i) - Tair(i))*(1/f)*UA_wall + InternalGains(t) + (Tsupply - Tair(i))*Cp_Air*Flow)*dt(t)/(n*AirCapacitance); %net change in temperature
        end
        Tsolid(i+1) = Tsolid(i) + ((Tdb(t) - Tsolid(i))*(1/(1-f))*UA_wall + (Tair(i+1) - Tsolid(i))*(1/f)*UA_wall + ExternalGains(t) + (Tsky^4 - (273+Tsolid(i))^4)*sigA).*dt(t)/(n*Build.VariableStruct.Capacitance*Build.Area);
    end
    Tzone(t+1) = Tair(end);
    Twall(t+1) = Tsolid(end);
end

end%Ends function building_simulate