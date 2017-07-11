function Out = ZonalBuilding(t,Y, Inlet,block,string1)
% a zonal HVAC building model, with air/liquid H2O heat exchangers for cooling and re-heating
% Eight (8) inlets: ambient air temperature, ambient air humidity, fan flow rate, HVAC temperature set point, dew point set point, damper position, cold water flow, hot water flow
% Five (5) outlets: Temperature, Humidity, Cooling Request, Heating Request, Mode (heating/cooling)
% Two*zones (2*n) states: Temperature, Humidity for each zone
global Tags
n = block.zones;
P = 101.325; %kPa, atmospheric pressure value
Cp_H2O = 4.184; % kJ/kg H2O

%ambient air
AmbientAir = makeAir(Inlet.Tamb,Inlet.Humidity,(1-Inlet.Damper)*Inlet.Flow,'abs');
%recirculated air
RecircAir = makeAir(Y(1:n),Y(n+1:end),Inlet.Damper*Inlet.Flow,'abs');
%mixed air
MixedAir = MixAir(RecircAir,AmbientAir);
Hmixed = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18);

%calculate treated air from cold and hot water flow 
TreatedAir = MixedAir;
Qcool = Inlet.CWflow.H2O*18*Cp_H2O.*((MixedAir.T - 2) - Inlet.CWflow.T); %assume cold water is heated up to within 2C of mixed air temp
P_H2O_mix = MixedAir.H2O/NetFlow(MixedAir)*P; %partial pressure of water wapor
mixed_DP = 6.54 + 14.526*log(P_H2O_mix) + 0.7389*log(P_H2O_mix)^2 + 0.09486*log(P_H2O_mix)^3 + 0.4569*(P_H2O_mix)^0.1984; %Dew point from partial pressure of water using ASHRAE 2013 Fundamentals eqn 39 valid from 0C to 93C

%find the temperature/H2O flow after cooling
error = ones(n,1);
while max(abs(error))>1e-3
    Hsupply = Inlet.Flow.*enthalpyAir(TreatedAir);%mass flow of dry air* enthalpy of dry air
    error = ((Hmixed-Qcool) - Hsupply)./(1.5*Inlet.Flow);
    TreatedAir.T = TreatedAir.T + error;
    for i = 1:1:block.zones
        if TreatedAir.T(i)<Inlet.CWflow.T(i)+2%check that air is not cooled to within 2C of the CW
            TreatedAir.T(i) = Inlet.CWflow.T(i)+2;
            if TreatedAir.T(i)<mixed_DP(i) %dehumidify if necessary
                dried = makeAir(TreatedAir.T(i),100,Inlet.Flow(i),'rel');
                TreatedAir.H2O(i) = dried.H2O;
            else
                TreatedAir.H2O(i) = MixedAir.H2O(i);
                dried.T = TreatedAir.T(i); dried.H2O = TreatedAir.H2O(i);dried.N2 = TreatedAir.N2(i);dried.O2 = TreatedAir.O2(i);dried.AR = TreatedAir.AR(i);dried.CO2 = TreatedAir.CO2(i);
            end
            Qcool(i) = Hmixed(i) - Inlet.Flow(i)*enthalpyAir(dried);
        elseif TreatedAir.T(i)<mixed_DP(i) %dehumidify if necessary
            dried = makeAir(TreatedAir.T(i),100,Inlet.Flow(i),'rel');
            TreatedAir.H2O(i) = dried.H2O;
        else TreatedAir.H2O(i) = MixedAir.H2O(i);
        end
    end
end
%find the enthalpy supplied by the hot water
Qheat = Inlet.HWflow.H2O*18*Cp_H2O.*(Inlet.HWflow.T - (TreatedAir.T + 2)); %assume hot water is cooled up to within 2C of treated air temp
    
if strcmp(string1,'Outlet')
    Out.Temperature = Y(1:n);
    Out.Humidity = Y(n+1:end);
    Out.Mode = Tags.(block.name).Mode;
    Tags.(block.name).Temperature = Y(1:n);
    Tags.(block.name).Humidity = Y(n+1:end);
    %calculate cooling and heating request
    sat_atDP = makeAir(Inlet.DPset,100,Inlet.Flow,'rel');
    DP_H2O = sat_atDP.H2O;
    
    %%cooling% cool to the supply temperature (if cooling is required)
    CooledAir = MixedAir;
    CooledAir.T = min(MixedAir.T,Inlet.Tset);%cooled mixed air
    %lower temperature and humidity for zones that need de-humidifying
    mixed_AH = MixedAir.H2O*18./(MassFlow(MixedAir)-MixedAir.H2O*18);%actual humidity: kg H20/kg dry air
    DHumid = (mixed_AH>block.maxHumidity);
    CooledAir.T(DHumid) = Inlet.DPset; 
    CooledAir.H2O(DHumid) = DP_H2O;
    %calculate cooling request
    Out.Qcool  = Hmixed - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    %calculate heating request
    HeatedAir = CooledAir;
    HeatedAir.T = Inlet.Tset;
    Out.Qheat  = enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    
    Out.CWreturn = Inlet.CWflow;
    Out.CWreturn.T = Inlet.CWflow.T + Qcool/(Inlet.CWflow.H2O*18*Cp_H2O);
    Out.HWreturn = Inlet.HWflow;
    Out.HWreturn.T = Inlet.HWflow.T + Qheat/(Inlet.HWflow.H2O*18*Cp_H2O);
elseif strcmp(string1,'dY')
    h = mod(t/3600,24);
    rho = interp1(linspace(35,-25,13),[1.1463 1.1653 1.1849 1.2052 1.2262 1.2479 1.2704 1.2938 1.3180 1.3432 1.3693 1.3965 1.4248],Y(1:n));%Density of dry air (kg/m^3)
    Occupancy = zeros(1,n);
    IntGains = zeros(1,n);
    for i = 1:1:n
        %% Occupancy and internal gains can be replaced with lookup functions of the smae name where the schedules are stored in SimSettings
        Occupancy(1,i) = interp1(linspace(0,24,length(block.OccupancySchedule(i,:))+1),block.Occupancy(i)*[block.OccupancySchedule(i,end),block.OccupancySchedule(i,:)]*block.Area(i),h);
        IntGains(1,i) = Occupancy(1,i)*120; %heat from occupants (W)
        IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.PlugSchedule(i,:))+1),block.PlugLoad(i)*[block.PlugSchedule(i,end),block.PlugSchedule(i,:)]*block.Area(i),h);%heat from plug loads (W)
        IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.LightingSchedule(i,:))+1),block.LightingLoad(i)*[block.LightingSchedule(i,end),block.LightingSchedule(i,:)]*block.Area(i),h);% Heat from lighting (W)
    end
    OutFlow = makeAir(Y(1:n),Y(n+1:end),Inlet.Flow,'abs');
    Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O).*enthalpyAir(OutFlow);
    Hin = Hsupply + Qheat;

    %calculate change in temperature & humidity states
    dY = 0*Y;
    dY(1:n) =(((Inlet.Tamb - Y(1:n))./block.Resistance) + IntGains + Hin*1000 - Hout*1000)./block.Capacitance;%Change in temperature
    dY(n+1:end) = (TreatedAir.H2O*18 + Occupancy*3e-6 - OutFlow.H2O*18)./(block.Volume.*rho);%Change in relative humidity = change in water mass/total air mass
    cool = (((Inlet.Tamb - Y(1:n))./block.Resistance + IntGains)>0);
    Tags.(block.name).Mode(cool) = -1;
    Tags.(block.name).Mode(~cool) = 1;
    Tags.(block.name).Cooling = Qcool; %kW of cooling provided
    Tags.(block.name).Heating = Qheat; %kW of heating provided
    Out = dY;
end