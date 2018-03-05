function [T_Heat,T_Cool] = EquilibTeperature(Build,Date,Weather,Twall,InternalGains,ExternalGains,Tmin,Tmax)
rho_Air = 1.225; % kg/m^3
n = 20;
nS = length(Date);
T_Heat = zeros(nS,1);
T_Cool = zeros(nS,1);
if nS == 1
    dt = 1;
else
    dt1 = Date(2) - Date(1);
    dt = (24*3600) * (Date - [Date(1)-dt1;Date(1:end-1)]); % duration of each time segment [seconds]
end
%% Ambient dewpoint
P = 101.325; % atmospheric pressure (kPa)
Tdb_K = Weather.Tdb+273.15; %Tdb (Kelvin)
satP = exp((-5.8002206e3)./Tdb_K + 1.3914993 - 4.8640239e-2*Tdb_K + 4.1764768e-5*Tdb_K.^2 - 1.4452093e-8*Tdb_K.^3 + 6.5459673*log(Tdb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
P_H2O = Weather.RH/100.*satP; % kPa
Cp_amb = 1.006 + 1.86*(.621945*(P_H2O./(P-P_H2O))); % kJ/kg*K
MinFlow = Build.VariableStruct.Volume*rho_Air*Build.VariableStruct.AirChangePerHr*(1/3600); %kg/s of flow
for t = 1:1:nS
    Trange = linspace(Tmin(t)-2*Build.VariableStruct.Comfort,Tmax(t)+2*Build.VariableStruct.Comfort,n)';
    NetGain = HVAC_equilibrium(Build,Twall(t),Weather.Tdb(t),Trange,InternalGains(t),ExternalGains(t),dt(t));
    needHeat = (min(Trange,Build.VariableStruct.ColdAirSet) - Weather.Tdb(t))*Cp_amb(t)*MinFlow*dt(t) - NetGain;
    if all(needHeat)>0
        T_Heat(t) = Trange(1);
    elseif all(needHeat)<0
        T_Heat(t) = Trange(end);
    else
        T_Heat(t) = interp1(Trange,needHeat,0);
    end
    needCool = (Weather.Tdb(t) - Trange)*Cp_amb(t)*MinFlow*dt(t) + NetGain;
    if all(needCool)>0
        T_Cool(t) = Trange(end);
    elseif all(needCool)<0
        T_Cool(t) = Trange(1);
    else
        T_Cool(t) = interp1(Trange,needCool,0);
    end
end
end%Ends function EquilibTemperature