function Out = HVACcyclic(t,Y, Inlet,block,string1)
% Controls for zonal building HVAC
% Four (4) inlets: temperature, cooling request, heating request, mode (heating or cooling)
% Six (6) outlets: Air mass flow, HVAC temperature set point, HVAC dew point, Damper position, cold water flow, hot water flow 
% Four (4) states: mass flow (of dry air), damper position, cooling valve position (cold water mass flow), heating valve position (hot water mass flow)
global Tags
n = block.zones;
Cp_H2O = 4.184; % kJ/kg H2O
Tol = block.Tolerance;
Terror = Inlet.Temperature - block.Target(1); %positive if room is too hot, negative if too cold
if strcmp(string1,'Outlet')
    Out.AirFlow = max(min(Y(1:n)+block.PropGain(1)*Terror,block.maxFlow),block.minFlow); 
    Out.Tset = block.ColdAirSetpoint*ones(block.zones,1);
    Out.Tset(Inlet.Mode == 1) = block.HotAirSetpoint;
    Out.DPset = block.DPset.IC; %unchanged
    Out.Damper = block.Damper.IC; 
    %%cold/hot water valve
    Out.ColdWater.T = block.ColdWater.IC.T;
    Out.HotWater.T = block.HotWater.IC.T;
    if Inlet.Mode == -1
        if Terror<-Tol ||  (Tags.(block.name).Cooling<1e-6 && Terror<Tol) %don't run AC
            Out.ColdWater.H2O = 0;
        else %if Terror>2 || (Tags.(block.name).Cooling>0 && Terror>-2)
            Out.ColdWater.H2O = block.CoolingPower./(Cp_H2O*18*(12)); %assumes 12 C of cooling, find kmol/s of H2O
        end
        Out.HotWater.H2O = Inlet.Qheat./(Cp_H2O*18*(Out.HotWater.T - Inlet.Temperature));%reheat
    else
        Out.ColdWater.H2O = 0;
        if Terror>Tol ||  (Tags.(block.name).Heating<1e-6 && Terror>-Tol) %don't run heater
            Out.HotWater.H2O = 0;
        else %if Terror<-2 || (Tags.(block.name).Heating>0 && Terror<2)
            Out.HotWater.H2O = block.HeatingPower./(Cp_H2O*18*(24)); %assumes 24 C of heating
        end
    end
    Tags.(block.name).FlowRate = Out.AirFlow;
    Tags.(block.name).Temperature = Out.Tset;
    Tags.(block.name).Cooling = Out.ColdWater.H2O*18*Cp_H2O*(12); %assumes 12 C of cooling
    Tags.(block.name).Heating = Out.HotWater.H2O*18*Cp_H2O*(24); %assumes 24 C of heating
    Tags.(block.name).CoolingPower = Tags.(block.name).Cooling/block.COP;
elseif strcmp(string1,'dY')
    %%mass flow%room is too cold, reduce flow rate, room is too hot, increase flow
    dY = block.Gain(1)*Terror;
    dY(Inlet.Mode==1) = (block.minFlow - Y(Inlet.Mode==1))/60;%heating mode (minimum flow)
    %anti-wind-up
    sat1 = nonzeros((1:n)'.*(dY>0 & Y>=1));
    dY(sat1) = (block.maxFlow - Y(sat1))/60;
    sat2 = nonzeros((1:n)'.*(dY<0 & Y<=block.minFlow/block.maxFlow));
    dY(sat2) = (block.minFlow - Y(sat2))/60;
    
    Out = dY;
end