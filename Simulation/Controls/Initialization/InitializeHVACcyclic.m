function block = InitializeHVACcyclic(varargin)
% Controls for zonal building HVAC with fixed power chiller / heater (reheat is variable)
% Four (4) inlets: temperature, cooling request, heating request, mode (heating or cooling)
% Six (6) outlets: Air mass flow, HVAC temperature set point, HVAC dew point, Damper position, cold water flow, hot water flow 
% Four (4) states: mass flow (of dry air), damper position, cooling valve position (cold water mass flow), heating valve position (hot water mass flow)
global Tags
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Mass flow (of dry air)';};
    
    block.InletPorts = {'Temperature','Qcool','Qheat','Mode'};
    block.Temperature.IC = 22.2; 
    block.Qcool.IC = 0; 
    block.Qheat.IC = 0; 
    block.Mode.IC = -1;

    block.OutletPorts = {'AirFlow';'Tset';'DPset';'Damper';'ColdWater';'HotWater';};
    block.AirFlow.IC = block.minFlow;
    block.Tset.IC = block.ColdAirSetpoint;
    block.DPset.IC = 11;
    block.Damper.IC = block.Entrainment/100;
    block.ColdWater.IC.T = 4;
    block.ColdWater.IC.H2O = 0;
    block.HotWater.IC.T = 40;
    block.HotWater.IC.H2O = 0;
    
    block.InitializeError = 1;
    
    block.Scale = [block.maxFlow;];
    block.IC = [block.minFlow/block.maxFlow;]; % inital condition for flow rate
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
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
    
    block.IC = flow/block.maxFlow;  % inital condition for flow rate, damper position (fraction of recirculated air), hot water flow, cold water flow
    block.InitializeError = abs(Terror); 
    Tags.(block.name).Temperature = block.Tset.IC;
    Tags.(block.name).Cooling = Inlet.Qcool;
    Tags.(block.name).Heating = Inlet.Qheat;
    Tags.(block.name).CoolingPower = Inlet.Qcool/block.COP;
end