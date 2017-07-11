function block = InitializeZonalBuilding(varargin)
% a zonal HVAC building model
% Eight (8) inlets: ambient air temperature, ambient air humidity, fan flow rate, HVAC temperature set point, dew point set point, damper position, cold water flow, hot water flow
% Five (5) outlets: Temperature, Humidity, Cooling Request, Heating Request, Mode (heating/cooling)
% Two*zones (2*n) states: Temperature, Humidity for each zone
global Tags
block = varargin{1};
if length(varargin)==1 % first initialization
    block.zones = length(block.Capacitance);
    block.Scale = [22.2*ones(block.zones,1); .0085*ones(block.zones,1)]; %temperature and humidity
    block.IC = ones(length(block.Scale),1);
    %%
    block.InletPorts = {'Tamb';'Humidity';'Flow';'Tset';'DPset';'Damper';'CWflow';'HWflow';};
    block.Tamb.IC = 25;
    block.Humidity.IC = 0.01;
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
    
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    t = 0;
    n = block.zones;
    Cp_H2O = 4.184; % kJ/kg H2O
    h = mod(t/3600,24);
    Occupancy = zeros(1,n);
    IntGains = zeros(1,n);
    for i = 1:1:n
        %% Occupancy and internal gains can be replaced with lookup functions of the smae name where the schedules are stored in SimSettings
        Occupancy(1,i) = interp1(linspace(0,24,length(block.OccupancySchedule(i,:))+1),block.Occupancy(i)*[block.OccupancySchedule(i,end),block.OccupancySchedule(i,:)]*block.Area(i),h);
        IntGains(1,i) = Occupancy(1,i)*120; %heat from occupants (W)
        IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.PlugSchedule(i,:))+1),block.PlugLoad(i)*[block.PlugSchedule(i,end),block.PlugSchedule(i,:)]*block.Area(i),h);%heat from plug loads (W)
        IntGains(1,i)  = IntGains(1,i)  + interp1(linspace(0,24,length(block.LightingSchedule(i,:))+1),block.LightingLoad(i)*[block.LightingSchedule(i,end),block.LightingSchedule(i,:)]*block.Area(i),h);% Heat from lighting (W)
    end
    %ambient air
    AmbientAir = makeAir(Inlet.Tamb,Inlet.Humidity,(1-Inlet.Damper)*Inlet.Flow,'abs');
    %recirculated air
    RecircAir = makeAir(block.Scale(1:n),block.Scale(n+1:end),Inlet.Damper*Inlet.Flow,'abs');
    %mixed air
    MixedAir = MixAir(RecircAir,AmbientAir);
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
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
    Tags.(block.name).Cooling = block.Qcool.IC; %kW of cooling provided
    Tags.(block.name).Heating = block.Qheat.IC; %kW of heating provided
    Tags.(block.name).Mode = block.Mode.IC;
end