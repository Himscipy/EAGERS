function block = InitializeHVAC_NREL(varargin)
% Controls for simple building HVAC temperaure control via mass flow
% Three (3) inlets: Building temperature, building mode (heating or cooling), ambient temperature
% Three (3) outlets: Supply flow to the building, supply temperature to the building Cooling load, heating load
% One (1) state: mass flow 
global Tags
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Air Supply to Building';};
    
    block.ColdAirSetpoint = 12.8; %minimum temperature of cooling air supplied by HVAC
    block.HotAirSetpoint = 40; %temperature of heating air supplied by HVAC when furnace is on
    block.Cp = 1e3; % J/kg*K
    
    block.InletPorts = {'Temperature','Mode','Tamb'};
    block.Temperature.IC = 22.2; 
    block.Mode.IC = 'cooling';
    block.Tamb.IC = 25;
        
    block.OutletPorts = {'massflow';'temperature';'Cooling';'Heating';};   
    block.massflow.IC = 1;
    block.temperature.IC = block.ColdAirSetpoint;
    block.Cooling.IC = 0;
    block.Heating.IC = 0;
    
    block.InitializeError = 1;
    
    block.Scale = 1;
    block.IC = 1; % inital condition 
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    %find temperature and flow
    Terror = Inlet.Temperature - block.Target(1);
    flow = max(1e-3,block.Scale(1));
    flow = max(0,(1 + Terror/10)*flow);
    Tmix = block.Entrainment/100*Inlet.Temperature + (1-block.Entrainment/100)*Inlet.Tamb;
    if strcmp(Inlet.Mode,'cooling')
        T = block.ColdAirSetpoint;
        block.Cooling.IC = flow*block.Cp*(Tmix - T);
        block.Heating.IC = 0;
    else %heating
        T = block.HotAirSetpoint;
        block.Cooling.IC = 0;
        block.Heating.IC = flow*block.Cp*(T - Tmix);
    end
   
    block.Scale = flow;
    block.massflow.IC = flow;
    block.temperature.IC = T;
    
    block.InitializeError = abs(Terror); %0.5*block.InitializeError + 0.5*
    Tags.(block.name).Cooling = block.Cooling.IC;
    Tags.(block.name).Heating = block.Heating.IC;
    Tags.(block.name).CoolingPower = block.Cooling.IC/block.COP;
end