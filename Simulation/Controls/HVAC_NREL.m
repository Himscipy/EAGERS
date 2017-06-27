function Out = HVAC_NREL(t,Y, Inlet,block,string1)
% Controls for simple building HVAC temperaure control via mass flow
% Three (3) inlets: Building temperature, building mode (heating or cooling), ambient temperature
% Three (3) outlets: Supply flow to the building, supply temperature to the building Cooling load, heating load
% One (1) state: mass flow 
global Tags
Terror = Inlet.Temperature - block.Target(1);
if strcmp(string1,'Outlet')
    Out.massflow = Y(1) + block.PropGain(1)*Terror;
    Tmix = block.Entrainment/100*Inlet.Temperature + (1-block.Entrainment/100)*Inlet.Tamb;
    if strcmp(Inlet.Mode,'cooling')
        Out.temperature = block.ColdAirSetpoint;
        Tags.(block.name).Cooling = Out.massflow*block.Cp*(Tmix - Out.temperature);
        Tags.(block.name).Heating = 0;
    else
        Out.temperature = block.HotAirSetpoint;
        Tags.(block.name).Cooling = 0;
        Tags.(block.name).Heating = Out.massflow*block.Cp*(Out.temperature - Tmix);
    end
    Tags.(block.name).FlowRate = Out.massflow;
    Tags.(block.name).Temperature = Out.temperature;
elseif strcmp(string1,'dY')
    dY = Y*0;
    %%mass flow
    if strcmp(Inlet.Mode,'cooling') %room is too cold, reduce flow rate, room is too hot, increase flow
        dY(1) = block.Gain(1)*Terror;
    else 
        dY(1) = -block.Gain(1)*Terror;
    end
    Out = dY;
end