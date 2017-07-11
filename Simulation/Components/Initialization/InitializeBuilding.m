function block = InitializeBuilding(varargin)
% a simple building model
% Five (5) inlets: air flow from HVAC system, treated air temperature (T after cooling, T after heating), damper position, ambient temperature, ambient humidity
% Three (3) outlets: Temperature, Humidity, Mode (heating/cooling)
% Two (2) state: Temperature, Humidity
global Tags
block = varargin{1};
if length(varargin)==1 % first initialization
    block.Scale = [22.2 .0085]; %temperature and humidity
    block.IC = ones(length(block.Scale),1);
    %%
    block.InletPorts = {'Flow','Tset','Damper','Tamb','ambHumidity'};
    block.Flow.IC = 1;
    block.Tset.IC = {12.8,12.8};
    block.Damper.IC = .75;
    block.Tamb.IC = 25;
    block.ambHumidity.IC = 0.01;
    
    block.OutletPorts = {'Temperature';'Humidity';'Mode';};
    block.Temperature.IC = 22.2;
    block.Humidity.IC = 0.0085;
    block.Mode.IC = 'cooling';
    
    block.P_Difference = {};
    
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
    Tags.(block.name).Mode = 'cooling';
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    P = 101.325; %ambient pressure in kPa
    t = 0;
    h = mod(t/3600,24);
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)
    
    %ambient air
    AmbientAir = makeAir(Inlet.Tamb,Inlet.ambHumidity,(1-Inlet.Damper)*Inlet.Flow,'abs');
    %recirculated air
    RecircAir = makeAir(block.Scale(1),block.Scale(2),Inlet.Damper*Inlet.Flow,'abs');
    %mixed air
    MixedAir = MixAir(RecircAir,AmbientAir);

    %%cooling
    CooledAir = MixedAir;
    if ~isempty(Inlet.Tset{1})
        CooledAir.T = Inlet.Tset{1}+273.15;% cool to the supply temperature if cooling 
        %check if it is below DP, if so remove moisture
        P_H2O = CooledAir.H2O/NetFlow(CooledAir)*P;
        satP = exp((-5.8002206e3)./CooledAir.T + 1.3914993 - 4.8640239e-2*CooledAir.T + 4.1764768e-5*CooledAir.T.^2 - 1.4452093e-8*CooledAir.T.^3 + 6.5459673*log(CooledAir.T))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
        if P_H2O>satP
            CooledAir.H2O = CooledAir.H2O * satP/P_H2O;
        end
    end
    Tags.(block.name).Cooling  = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    %%heating
    HeatedAir = CooledAir;
    if ~isempty(Inlet.Tset{2})
        HeatedAir.T = Inlet.Tset{2}+273.15;
    end
    Tags.(block.name).Heating  = enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);

    Hin = (MassFlow(HeatedAir) - 18*HeatedAir.H2O)*enthalpyAir(HeatedAir); %mass flow of dry air* enthalpy of dry air
    OutFlow = HeatedAir;
    OutFlow.H2O = OutFlow.H2O + Occupancy*3e-6;%3x10^-6 kg/s is the average rate of water vapor released per person due to respiration  

    error = 1;
    while abs(error)>1e-3
        Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O)*enthalpyAir(OutFlow);
        Herror = 1/(1000*block.Resistance)*(Inlet.Tamb + block.Resistance*(IntGains + Hin*1000 - Hout*1000) + 273.15 - OutFlow.T); %temperature which  balances energy equation
        OutFlow.T = OutFlow.T + Herror/(1.5*MassFlow(OutFlow));
        error = (block.Scale(1) - (OutFlow.T-273.15));
        block.Scale(1) = OutFlow.T-273.15; %temperature state is in celcius
    end
    
    block.Scale(2) = OutFlow.H2O*18/(MassFlow(OutFlow)-OutFlow.H2O*18);%absolute humidity
    block.Temperature.IC = block.Scale(1);
    block.Humidity.IC = block.Scale(2);
    if ((Inlet.Tamb - (OutFlow.T-273.15))/block.Resistance + IntGains)>0
        block.Mode.IC = 'cooling';
    else
        block.Mode.IC = 'heating';
    end
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).Humidity = block.Humidity.IC;
    Tags.(block.name).Mode = block.Mode.IC;
end