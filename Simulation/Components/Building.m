function Out = Building(t,Y, Inlet,block,string1)
% a simple building model
% Five (5) inlets: air flow from HVAC system, treated air temperature (T after cooling, T after heating), damper position, ambient temperature, ambient humidity
% Three (3) outlets: Temperature, Humidity, Mode (heating/cooling)
% Two (2) state: Temperature, Humidity
global Tags    
if strcmp(string1,'Outlet')
    Out.Humidity = Y(2);
    Out.Temperature = Y(1);
    
    Out.Mode = Tags.(block.name).Mode;
    Tags.(block.name).Temperature = Y(1);
    Tags.(block.name).Humidity = Y(2);
elseif strcmp(string1,'dY')
    h = mod(t/3600,24);
    rho = interp1(linspace(35,-25,13),[1.1463 1.1653 1.1849 1.2052 1.2262 1.2479 1.2704 1.2938 1.3180 1.3432 1.3693 1.3965 1.4248],Y(1));%Density of dry air (kg/m^3)
    P = 101.325; %ambient pressure
    %ambient air
    AmbientAir = makeAir(Inlet.Tamb,Inlet.ambHumidity,(1-Inlet.Damper)*Inlet.Flow,'abs');
    %recirculated air
    RecircAir = makeAir(Y(1),Y(2),Inlet.Damper*Inlet.Flow,'abs');
    %mixed air
    MixedAir = MixAir(RecircAir,AmbientAir);

    %%cooling
    CooledAir = MixedAir;
    if ~isempty(Inlet.Tset{1})
        CooledAir.T = Inlet.Tset{1}+273.15;% cool to the supply temperature if cooling 
        %check if it is below DP, if so remove moisture
        P_H2O = CooledAir.H2O/NetFlow(CooledAir);
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

    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)

    OutFlow = makeAir(Y(1),Y(2),Inlet.Flow,'abs');%dry air mass flow leaving is equal to the dry air mass flow entering
    
    Hin = (MassFlow(HeatedAir) - 18*HeatedAir.H2O)*enthalpyAir(HeatedAir); %mass flow of dry air* enthalpy of dry air
    Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O)*enthalpyAir(OutFlow);
    
    dY = 0*Y;
    dY(1) =(((Inlet.Tamb - Y(1))/block.Resistance) + IntGains + Hin*1000 - Hout*1000)/block.Capacitance;%Change in temperature
    dY(2) = (HeatedAir.H2O*18 + Occupancy*3e-6 - OutFlow.H2O*18)/(block.Volume*rho);%Change in relative humidity = change in water mass/total air mass
    if ((Inlet.Tamb - Y(1))/block.Resistance + IntGains)>0
        Tags.(block.name).Mode = 'cooling';
    else Tags.(block.name).Mode = 'heating';
    end
    Out = dY;
end