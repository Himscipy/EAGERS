function Out = ZonalBuilding(varargin)
% a zonal HVAC building model, with air/liquid H2O heat exchangers for cooling and re-heating
% Eight (8) inlets: ambient air temperature, ambient air humidity, fan flow rate, HVAC temperature set point, dew point set point, damper position, cold water flow, hot water flow
% Five (5) outlets: Temperature, Humidity, Cooling Request, Heating Request, Mode (heating/cooling)
% Two*zones (2*n) states: Temperature, Humidity for each zone
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.zones = length(block.Capacitance);
    block.Scale = [22.2*ones(block.zones,1); .0085*ones(block.zones,1)]; %temperature and humidity
    block.IC = ones(length(block.Scale),1);
    block.UpperBound = inf*ones(length(block.Scale),1);
    block.LowerBound = zeros(length(block.Scale),1);
    %%
    block.InletPorts = {'Tamb';'ambHumidity';'Flow';'Tset';'DPset';'Damper';'CWflow';'HWflow';};
    block.Tamb.IC = 25;
    block.ambHumidity.IC = 0.01;
    block.Flow.IC = 1*ones(block.zones,1); %mass flow of dry air supply
    block.Tset.IC = 12.8*ones(block.zones,1); %zone HVAC temperature set point
    block.DPset.IC = 11*ones(block.zones,1); %dew point set point
    block.Damper.IC = 0.4*ones(block.zones,1); %split between recirculated air and fresh air
    block.CWflow.IC.T = 4*ones(block.zones,1); %cold water temperature
    block.CWflow.IC.H2O = 0*ones(block.zones,1); %cold water flow
    block.HWflow.IC.H2O = 40*ones(block.zones,1); %hot water temperature
    block.HWflow.IC.H2O = 0*ones(block.zones,1); %hot water flow
    
    block.OutletPorts = {'Temperature';'Humidity';'Qcool';'Qheat';'CWreturn';'HWreturn';'Mode';};
    block.Temperature.IC = 22.2*ones(block.zones,1); %temperature of each zone
    block.Humidity.IC = 0.0085*ones(block.zones,1); %absolute humidity of each zone kg H2O/ kg dry air
    block.Qcool.IC = 0*ones(block.zones,1);%request for cooling (kW)
    block.Qheat.IC = 0*ones(block.zones,1);%request for heating (kW)
    block.CWreturn.IC.T = 18*ones(block.zones,1); %cold water return temperature
    block.CWreturn.IC.H2O = 0*ones(block.zones,1); %cold water return flow
    block.HWreturn.IC.H2O = 20*ones(block.zones,1); %hot water return temperature
    block.HWreturn.IC.H2O = 0*ones(block.zones,1); %hot water return flow
    block.Mode.IC = -1*ones(block.zones,1); %-1 = cooling, 1 = heating
    
    block.P_Difference = {};
    
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    t = 0;
    n = block.zones;
    Cp_H2O = 4.184; % kJ/kg H2O
    [MixedAir, IntGains, Occupancy] = findGains_andMixedAir(Inlet,block,t,block.Scale);
    mixed_AH = MixedAir.H2O*18./(MassFlow(MixedAir)-MixedAir.H2O*18);%actual humidity: kg H20/kg dry air
    sat_atDP = makeAir(Inlet.DPset,100,Inlet.Flow,'rel');
    DP_H2O = sat_atDP.H2O;
    
    %%cooling% cool to the supply temperature (if cooling is required)
    CooledAir = MixedAir;
    CooledAir.T = min(MixedAir.T,Inlet.Tset+273.15);%cooled mixed air
    %lower temperature and humidity for zones that need de-humidifying
    DHumid = (mixed_AH>block.maxHumidity);
    CooledAir.T(DHumid) = Inlet.DPset+273.15; 
    CooledAir.H2O(DHumid) = DP_H2O;
    HeatedAir = CooledAir;
    HeatedAir.T = Inlet.Tset+273.15;
 
    %equilibrium humidity (HVAC humidity + Occupant H2O flow rate) / dry air flow rate
    block.Scale(n+1:end) = (HeatedAir.H2O*18 + Occupancy*3e-6)/Inlet.Flow; %3x10^-6 kg/s is the average rate of water vapor released per person due to respiration  
    OutFlow = makeAir(block.Scale(1:n),block.Scale(n+1:end),Inlet.Flow,'abs');
    
    %during initialization we can assume we recieved the heating/cooling request to hit the HVAC set point
    Hin = (MassFlow(HeatedAir) - 18*HeatedAir.H2O).*enthalpyAir(HeatedAir); %mass flow of dry air* enthalpy of dry air
    error = 1;
    while max(abs(error))>1e-3
        Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O).*enthalpyAir(OutFlow);
        OutFlow.T = (- Inlet.Tamb + block.Resistance.*(IntGains + Hin*1000 - Hout*1000)) + 273.15; %temperature which  balances energy equation
        error = (block.Scale(1:n) - (OutFlow.T-273.15));
        block.Scale(1:n) = OutFlow.T-273.15; %temperature state is in celcius
    end
    
    block.Temperature.IC = block.Scale(1:n);
    block.Humidity.IC = block.Scale(n+1:end);
    %calculate cooling request
    block.Qcool.IC  = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    %calculate heating request
    block.Qheat.IC  = enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    %calculate return water
    block.CWreturn.IC = Inlet.CWflow;
    block.CWreturn.IC.T = Inlet.CWflow.T + block.Qcool.IC/(Inlet.CWflow.H2O*18*Cp_H2O);
    block.HWreturn.IC = Inlet.HWflow;
    block.HWreturn.IC.T = Inlet.HWflow.T + block.Qheat.IC/(Inlet.HWflow.H2O*18*Cp_H2O);
    cool = (((Inlet.Tamb - (OutFlow.T-273.15))/block.Resistance + IntGains)>0);
    block.Mode.IC(cool) = -1;
    block.Mode.IC(~cool) = 1;
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    n = block.zones;
    P = 101.325; %kPa, atmospheric pressure value
    Cp_H2O = 4.184; % kJ/kg H2O
    [MixedAir, IntGains, Occupancy] = findGains_andMixedAir(Inlet,block,t,Y);
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
    cool = (((Inlet.Tamb - Y(1:n))./block.Resistance + IntGains)>0);
    if strcmp(string1,'Outlet')
        Out.Temperature = Y(1:n);
        Out.Humidity = Y(n+1:end);
        Out.Mode(cool) = -1;
        Out.Mode(~cool) = 1;
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
        Tags.(block.name).Temperature = Y(1:n);
        Tags.(block.name).Humidity = Y(n+1:end);
        Tags.(block.name).Mode = Out.Mode;
        Tags.(block.name).Cooling = Qcool; %kW of cooling provided
        Tags.(block.name).Heating = Qheat; %kW of heating provided
    elseif strcmp(string1,'dY')
        rho = interp1(linspace(35,-25,13),[1.1463 1.1653 1.1849 1.2052 1.2262 1.2479 1.2704 1.2938 1.3180 1.3432 1.3693 1.3965 1.4248],Y(1:n));%Density of dry air (kg/m^3)
        OutFlow = makeAir(Y(1:n),Y(n+1:end),Inlet.Flow,'abs');
        Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O).*enthalpyAir(OutFlow);
        Hin = Hsupply + Qheat;

        %calculate change in temperature & humidity states
        dY = 0*Y;
        dY(1:n) =(((Inlet.Tamb - Y(1:n))./block.Resistance) + IntGains + Hin*1000 - Hout*1000)./block.Capacitance;%Change in temperature
        dY(n+1:end) = (TreatedAir.H2O*18 + Occupancy*3e-6 - OutFlow.H2O*18)./(block.Volume.*rho);%Change in relative humidity = change in water mass/total air mass
        
        
        Out = dY;
    end
end
end%Ends function ZonalBuilding

function [MixedAir, IntGains, Occupancy] = findGains_andMixedAir(Inlet,block,t,Y)
h = mod(t/3600,24);
n = block.zones;
%ambient air
AmbientAir = makeAir(Inlet.Tamb,Inlet.ambHumidity,(1-Inlet.Damper)*Inlet.Flow,'abs');
%recirculated air
RecircAir = makeAir(Y(1:n),Y(n+1:end),Inlet.Damper*Inlet.Flow,'abs');
%mixed air
MixedAir = MixAir(RecircAir,AmbientAir);
Occupancy = zeros(1,n);
IntGains = zeros(1,n);
for i = 1:1:n
    %% Occupancy and internal gains can be replaced with lookup functions of the smae name where the schedules are stored in SimSettings
    Occupancy(1,i) = interp1(linspace(0,24,length(block.OccupancySchedule(i,:))+1),block.Occupancy(i)*[block.OccupancySchedule(i,end),block.OccupancySchedule(i,:)]*block.Area(i),h);
    IntGains(1,i) = Occupancy(1,i)*120; %heat from occupants (W)
    IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.PlugSchedule(i,:))+1),block.PlugLoad(i)*[block.PlugSchedule(i,end),block.PlugSchedule(i,:)]*block.Area(i),h);%heat from plug loads (W)
    IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.LightingSchedule(i,:))+1),block.LightingLoad(i)*[block.LightingSchedule(i,end),block.LightingSchedule(i,:)]*block.Area(i),h);% Heat from lighting (W)
end
end%Ends function findGains_andMixedAir