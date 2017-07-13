function Out = HVACcyclic(varargin)
% Controls for zonal building HVAC
% Four (4) inlets: temperature, cooling request, heating request, mode (heating or cooling)
% Six (6) outlets: Air mass flow, HVAC temperature set point, HVAC dew point, Damper position, cold water flow, hot water flow 
% Four (4) states: mass flow (of dry air), damper position, cooling valve position (cold water mass flow), heating valve position (hot water mass flow)
global Tags
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Mass flow (of dry air)';};
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
    block.InletPorts(end+1:end+4) = {'Temperature','Qcool','Qheat','Mode'};
    block.Temperature.IC = 22.2;
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
    
    block.Scale = [block.maxFlow;];
    block.IC = [block.minFlow/block.maxFlow;]; % inital condition for flow rate
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    
    %%mass flow & temperature
    Terror = Inlet.Temperature - block.Target(1); %posiive if room is too hot, negative if too cold
    flow = block.AirFlow.IC; %mass flow of dry air
    block.Tset.IC(Inlet.Mode==1) = block.HotAirSetpoint;
    sat1 = (Inlet.Mode==-1 & (Terror<0 & flow-block.minFlow<1e-10) | (Terror>0 & block.maxFlow-flow<1e-10));
    block.Tset.IC(sat1) = block.Tset.IC(sat1) - Terror(sat1);
    A = (~sat1 & Inlet.Mode==-1);
    block.Tset.IC(A) = block.ColdAirSetpoint;
    flow(A) = min(block.maxFlow,max(block.minFlow,(1 - Terror(A)/5).*flow(A)));
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
    
    block.IC = flow/block.maxFlow;  % inital condition for flow rate, damper position (fraction of recirculated air), hot water flow, cold water flow
    block.InitializeError = abs(Terror); 
    Tags.(block.name).Temperature = block.Tset.IC;
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
    Tol = block.Tolerance;
    Terror = Inlet.Temperature - block.Target(1); %positive if room is too hot, negative if too cold
    if strcmp(string1,'Outlet')
        Out.Measured1 = Inlet.Temperature;
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
end
end%Ends function HVACcyclic