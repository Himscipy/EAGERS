function Out = HVACproportional(t,Y, Inlet,block,string1)
% Controls for zonal building HVAC
% Four (4) inlets: temperature, cooling request, heating request, mode (heating or cooling)
% Six (6) outlets: Air mass flow, HVAC temperature set point, HVAC dew point, Damper position, cold water flow, hot water flow 
% Four (4) states: mass flow (of dry air), damper position, cooling valve position (cold water mass flow), heating valve position (hot water mass flow)
global Tags
n = block.zones;
Cp_H2O = 4.184; % kJ/kg H2O
% Tol = block.Tolerance;
Terror = Inlet.Temperature - block.Target(1); %positive if room is too hot, negative if too cold

Out.ColdWater.T = block.ColdWater.IC.T;
Out.HotWater.T = block.HotWater.IC.T;
ColdWater = Inlet.Qcool/(Cp_H2O*18*(Inlet.Temperature - Out.ColdWater.T));
CWerror = ColdWater - (Y(3) + block.PropGain(3)*ColdWater)/(1+block.PropGain(3));
HotWater = Inlet.Qheat/(Cp_H2O*18*(Out.HotWater.T - Inlet.Temperature));
HWerror = HotWater - (Y(4) + block.PropGain(4)*HotWater)/(1+block.PropGain(4));
if strcmp(string1,'Outlet')
    Out.AirFlow = max(min(Y(1:n)+block.PropGain(1)*Terror,block.maxFlow),block.minFlow); 
    Out.Tset = block.ColdAirSetpoint*ones(block.zones,1);
    Out.Tset(Inlet.Mode == 1) = block.HotAirSetpoint;
    Out.DPset = block.DPset.IC; %unchanged
    Out.Damper = Y(2); 
    %%cold/hot water valve
    Out.ColdWater.H2O = Y(3) + block.PropGain(3)*CWerror;
    Out.HotWater.H2O = Y(4) + block.PropGain(4)*HWerror;
    Tags.(block.name).FlowRate = Out.AirFlow;
    Tags.(block.name).Temperature = Out.Tset;
    Tags.(block.name).Cooling = Out.ColdWater.H2O*18*Cp_H2O*(12); %assumes 12 C of cooling
    Tags.(block.name).Heating = Out.HotWater.H2O*18*Cp_H2O*(24); %assumes 24 C of heating
    Tags.(block.name).CoolingPower = Tags.(block.name).Cooling/block.COP;
elseif strcmp(string1,'dY')
    %%mass flow%room is too cold, reduce flow rate, room is too hot, increase flow
    dY1 = block.Gain(1)*Terror;
    Y1 = Y(1:n);
    dY1(Inlet.Mode==1) = (block.minFlow - Y1(Inlet.Mode==1))/60;%heating mode (minimum flow)
    %anti-wind-up
    sat1 = nonzeros((1:n)'.*(dY1>0 & Y1>=1));
    dY1(sat1) = (block.maxFlow - Y1(sat1))/60;
    sat2 = nonzeros((1:n)'.*(dY1<0 & Y1<=block.minFlow/block.maxFlow));
    dY1(sat2) = (block.minFlow - Y1(sat2))/60;
    
    %%damper position
    MixedAir.T = Inlet.Tamb*Y(n+1:2*n) + Inlet.Temperature*(1-Y(n+1:2*n));%estimate of mixed air temperature
    dY2 = block.Gain(2)*(MixedAir.T - block.Target(1));
    if Inlet.Mode ==1
        dY2(Inlet.Mode ==1) = .25 - block.Gain(2)*Y(n+1:2*n); %heating, close damper
    end
    %anti-wind-up
    sat1 = nonzeros((1:n)'.*(dY2>0 & Y(n+1:2*n)>=1));
    dY2(sat1) = 1-block.Gain(2)*Y(sat1+n);
    sat2 = nonzeros((1:n)'.*(dY2<0 & Y(n+1:2*n)<=0.25));
    dY2(sat2) = 0.25-block.Gain(2)*Y(sat2+n);

    %%cold water valve
    dY3 = block.Gain(3)*(ColdWater-Y(2*n+1:3*n))/10; % /10 is the time constant (10 s)
    %%hot water valve
    dY4 = block.Gain(4)*(HotWater-Y(3*n+1:4*n))/10; % /10 is the time constant (10 s)
    Out = [dY1;dY2;dY3;dY4;];
end