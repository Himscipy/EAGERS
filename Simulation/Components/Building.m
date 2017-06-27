function Out = Building(t,Y, Inlet,block,string1)
% a simple building model
% Two (2) inlets: air flow from HVAC system,  ambient temperature
% Two (2) outlets: Temperature, Humidity
% Two (2) state: Temperature, Humidity
global Tags
if Y(1)>40 || Y(1)<15
%     disp('WTF')
%something went wrong, use previous to avoid error
    Y(1) = Tags.(block.name).Temperature;
    Y(2) = Tags.(block.name).Humidity;
end
if strcmp(string1,'Outlet')
    Out.Humidity = Y(2);
    Out.Temperature = Y(1);
    
    Out.Mode = Tags.(block.name).Mode;
    Tags.(block.name).Temperature = Y(1);
    Tags.(block.name).Humidity = Y(2);
elseif strcmp(string1,'dY')
    h = mod(t/3600,24);
    rho = interp1(linspace(35,-25,13),[1.1463 1.1653 1.1849 1.2052 1.2262 1.2479 1.2704 1.2938 1.3180 1.3432 1.3693 1.3965 1.4248],Y(1));%Density of dry air (kg/m^3)
    Occupancy = interp1(linspace(0,24,length(block.OccupancySchedule)+1),block.Occupancy*[block.OccupancySchedule(end),block.OccupancySchedule]*block.Area,h);

    %dry air mass flow leaving is equal to the dry air mass flow entering + the mass flow of water from the occupants
    flowOut = MassFlow(Inlet.HVAC) - Inlet.HVAC.H2O*18 + Occupancy*3e-6;%3x10^-6 kg/s is the average rate of water vapor released per person due to respiration
    OutFlow = makeAir(Y(1),Y(2),flowOut,'abs');
    
    Hin = (MassFlow(Inlet.HVAC) - 18*Inlet.HVAC.H2O)*enthalpyAir(Inlet.HVAC); %mass flow of dry air* enthalpy of dry air
    Hout = (MassFlow(OutFlow) - 18*OutFlow.H2O)*enthalpyAir(OutFlow);
    
    IntGains = Occupancy*120; %heat from occupants (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.PlugSchedule)+1),block.PlugLoad*[block.PlugSchedule(end),block.PlugSchedule]*block.Area,h);%heat from plug loads (W)
    IntGains  = IntGains  + interp1(linspace(0,24,length(block.LightingSchedule)+1),block.LightingLoad*[block.LightingSchedule(end),block.LightingSchedule]*block.Area,h);% Heat from lighting (W)

    dY = 0*Y;
    dY(1) =(((Inlet.Tamb - Y(1))/block.Resistance) + IntGains + Hin*1000 - Hout*1000)/block.Capacitance;%Change in temperature
    dY(2) = (Inlet.HVAC.H2O*18 + Occupancy*3e-6 - OutFlow.H2O*18)/(block.Volume*rho);%Change in relative humidity = change in water mass/total air mass
    if ((Inlet.Tamb - Y(1))/block.Resistance + IntGains)>0
        Tags.(block.name).Mode = 'cooling';
    else Tags.(block.name).Mode = 'heating';
    end
    Out = dY;
end