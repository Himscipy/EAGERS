function Out = HVACproportional(varargin)
% Controls for zonal building HVAC
% Four (4) inlets: temperature, cooling request, heating request, mode (heating or cooling)
% Six (6) outlets: Air mass flow, HVAC temperature set point, HVAC dew point, Damper position, cold water flow, hot water flow 
% Four (4) states: mass flow (of dry air), damper position, cooling valve position (cold water mass flow), heating valve position (hot water mass flow)
global Tags
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Mass flow (of dry air)';'Damper position';'Cooling valve position (cold water mass flow)';'Heating valve position (hot water mass flow)';};
    block.TargetDescription = {'Building Temperature Setpoint (drybulb)';'Dewpoint Setpoint';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(block.Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target = Target;
    
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+4) = {'Temperature','Tamb','Qcool','Qheat','Mode'};
    block.Temperature.IC = 22.2;
    block.Tamb.IC = 25;
    block.Qcool.IC = 0; 
    block.Qheat.IC = 0; 
    block.Mode.IC = -1;
    
    block.OutletPorts(end+1:end+6) = {'AirFlow';'Tset';'DPset';'Damper';'ColdWater';'HotWater';};
    block.AirFlow.IC = block.minFlow;
    block.Tset.IC = block.ColdAirSetpoint;
    block.DPset.IC = 11;
    block.Damper.IC = block.Entrainment/100;
    block.ColdWater.IC.T = 4;
    block.ColdWater.IC.H2O = 0;
    block.HotWater.IC.T = 40;
    block.HotWater.IC.H2O = 0;
    
    block.Scale = [block.maxFlow; 1;1;1;];
    block.IC = [block.minFlow/block.maxFlow; block.Damper.IC; block.ColdWater.IC.H2O; block.HotWater.IC.H2O;]; % inital condition for flow rate, damper position (fraction of recirculated air), heating coil valve, cooling coil valve
    block.UpperBound = [inf,1,1,1];
    block.LowerBound = [0,0,0,0];
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Terror = Inlet.Temperature - block.Target(1); %posiive if room is too hot, negative if too cold
    flow = block.AirFlow.IC; %mass flow of dry air
    
    %%mass flow & temperature
    block.Tset.IC(Inlet.Mode==1) = block.HotAirSetpoint;
    sat1 = (Inlet.Mode==-1 & (Terror<0 & flow-block.minFlow<1e-10) | (Terror>0 & block.maxFlow-flow<1e-10));
    block.Tset.IC(sat1) = block.Tset.IC(sat1) - Terror(sat1);
    A = (~sat1 & Inlet.Mode==-1);
    block.Tset.IC(A) = block.ColdAirSetpoint;
    flow(A) = min(block.maxFlow,max(block.minFlow,(1 - Terror(A)/(Inlet.Tamb-Inlet.Temperature)).*flow(A)));
    block.AirFlow.IC = flow;
    %%dew point & damper position (unchanged)
    
    %%cold/hot water valve
    Cp_H2O = 4.184; % kJ/kg H2O
    block.ColdWater.IC.H2O = Inlet.Qcool/(Cp_H2O*18*(Inlet.Temperature - block.ColdWater.IC.T));
    block.HotWater.IC.H2O = Inlet.Qheat/(Cp_H2O*18*(block.HotWater.IC.T - Inlet.Temperature));
    
    block.Measured1.IC = Inlet.Temperature;
%     P = 101.325;
%     P_H2O = P*Inlet.Humidity/(0.621945+Inlet.Humidity);
%     T_K = Inlet.Temperature;
%     satP = exp((-5.8002206e3)./T_K + 1.3914993 - 4.8640239e-2*T_K + 4.1764768e-5*T_K.^2 - 1.4452093e-8*T_K.^3 + 6.5459673*log(T_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
%     RH = P_H2O/satP*100;
%     T_WB = Inlet.Temperature*atan(0.151977*(RH + 8.313659)^.5) + atan(Inlet.Temperature + RH) - atan(RH - 1.676331) + 0.00391838*(RH)^(3/2)*atan(0.023101*RH) - 4.686035; % http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
%     block.Measured2.IC = T_WB;
    
    block.IC = [flow/block.maxFlow; block.Damper.IC; block.ColdWater.IC.H2O; block.HotWater.IC.H2O;]; % inital condition for flow rate, damper position (fraction of recirculated air), hot water flow, cold water flow
    block.InitializeError = abs(Terror); 
    Tags.(block.name).Cooling = Inlet.Qcool;
    Tags.(block.name).Heating = Inlet.Qheat;
    Tags.(block.name).CoolingPower = Inlet.Qcool/block.COP;
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
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
        Out.Measured1 = Inlet.Temperature;
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
end
end%Ends function HVACproportional