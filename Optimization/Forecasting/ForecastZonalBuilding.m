function [Cooling, Heating] = ForecastZonalBuilding(varargin)
% a zonal HVAC building forecast to calculate the cooling and heating
% necessary to achieve the desired temperature (Tset) and the times specified in Date
% build is a structure of building parameters
% IC is a structure of the current time and current temperature
% Date is a mx1 vector of times at which you are requesting the cooling & heating power
% Tset is a m x n matrix of temperature settings
% DPset is a m x n matrix of dew point settings
% Damper is a m x n matrix of damper position (ambient air entrainment) settings
% Tamb is a m x 1 vector of ambient temperature
% Tamb is a m x 1 vector of ambient humidity
build = varargin{1};
IC = varargin{2};
Date = varargin{3};
Tamb = varargin{4};
Tset = varargin{5};
if length(varargin)>5
    Flow =  varargin{6};
    DPset = varargin{7};
    Damper = varargin{8};
    Humidity = varargin{9};
end
n = build.zones;
t = 86400*(Date-floor(Date)); %time in seconds of a 24 hour day to lookup internal gains and occupancy
dt = 86400*(Date - [IC.Timestamp;Date(1:end-1)]); %duration of each time segment
h = mod(t/3600,24);
ColdAirSet = 12.8;
Occupancy = zeros(length(t),n);
IntGains = zeros(length(t),n);
Cooling = zeros(length(t),n);
Heating = zeros(length(t),n);

for i = 1:1:n
    %% Occupancy and internal gains can be replaced with lookup functions of the same name where the schedules are stored in SimSettings
    Occupancy(:,i) = interp1(linspace(0,24,length(build.OccupancySchedule(i,:))+1),build.Occupancy(i)*[build.OccupancySchedule(i,end),build.OccupancySchedule(i,:)]*build.Area(i),h);
    IntGains(:,i) = Occupancy(:,i)*120; %heat from occupants (W)
    IntGains(:,i)  = IntGains(:,i)  + interp1(linspace(0,24,length(build.PlugSchedule(i,:))+1),build.PlugLoad(i)*[build.PlugSchedule(i,end),build.PlugSchedule(i,:)]*build.Area(i),h);%heat from plug loads (W)
    IntGains(:,i)  = IntGains(:,i)  + interp1(linspace(0,24,length(build.LightingSchedule(i,:))+1),build.LightingLoad(i)*[build.LightingSchedule(i,end),build.LightingSchedule(i,:)]*build.Area(i),h);% Heat from lighting (W)
    
    Ej = ((Tamb - Tset(:,i))./build.Resistance + IntGains(:,i)).*dt + (Tset(:,i) - [IC.T(i);Tset(1:end-1,i)])*build.Capacitance; %energy requirement in J
    Reheat = zeros(length(t));
    if ~isempty(Humidity)
        %ambient air
        AmbientAir = makeAir(Tamb,Humidity,(1-Damper(:,i))*Flow,'abs');
        %recirculated air
        RecircAir = makeAir(Tset,DPset,Damper(:,i)*Flow,'dp');
        %mixed air
        MixedAir = MixAir(RecircAir,AmbientAir);
        mixed_AH = MixedAir.H2O*18./(MassFlow(MixedAir)-MixedAir.H2O*18);%actual humidity: kg H20/kg dry air

        sat_atDP = makeAir(DPset(:,i),100,1,'rel');
        DP_H2O = sat_atDP.H2O;
        
        MixedAir.T = min(MixedAir.T,ColdAirSet+273.15);%cooled mixed air, this is the standard cooling calculated as Ej
        
        %lower temperature and humidity for zones that need de-humidifying
        CooledAir = MixedAir;
        DHumid = (mixed_AH>build.maxHumidity);
        CooledAir.T(DHumid) = Inlet.DPset+273.15; 
        CooledAir.H2O(DHumid) = DP_H2O; %remove water
        %calculate reheat load
        cool  = (enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18))*1000.*dt; %dehumidification energy in J
        
        HeatedAir = CooledAir;
        HeatedAir.T = ColdAirSet+273.15;
        %calculate reheat load
        Reheat  = (enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18))*1000.*dt; %reheat energy in J
    end
    Cooling(Ej<0,i) = Ej(Ej<0) + cool(Ej<0);
    Heating(Ej>0,i) = Ej(Ej>0) + Reheat(Ej<0);
end
end%Ends function ForecastZonal Building   
