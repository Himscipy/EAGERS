function Out = HVAC_NREL(varargin)
% Controls for simple building HVAC temperaure control via mass flow
% Three (3) inlets: Building temperature, building mode (heating or cooling), ambient temperature
% Three (3) outlets: Supply flow to the building, supply temperature to the building Cooling load, heating load
% One (1) state: mass flow 
global Tags
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Mass flow (of dry air)';};
    block.TargetDescription = {'Building Temperature Setpoint (drybulb)';};
    
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
    block.InletPorts(end+1:end+3) = {'Temperature','Mode','Tamb'};
    block.Temperature.IC = 22.2; 
    block.Mode.IC = 'cooling';
    block.Tamb.IC = 25;
    
    block.OutletPorts(end+1:end+4) = {'massflow';'temperature';'Cooling';'Heating';};   
    block.massflow.IC = 1;
    block.temperature.IC = block.ColdAirSetpoint;
    block.Cooling.IC = 0;
    block.Heating.IC = 0;
    
    block.ColdAirSetpoint = 12.8; %minimum temperature of cooling air supplied by HVAC
    block.HotAirSetpoint = 40; %temperature of heating air supplied by HVAC when furnace is on
    block.Cp = 1e3; % J/kg*K

    block.Scale = 1;
    block.IC = 1; % inital condition 
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
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
    block.Measured1.IC = Inlet.Temperature;
    
    block.InitializeError = abs(Terror); %0.5*block.InitializeError + 0.5*
    Tags.(block.name).Cooling = block.Cooling.IC;
    Tags.(block.name).Heating = block.Heating.IC;
    Tags.(block.name).CoolingPower = block.Cooling.IC/block.COP;
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Terror = Inlet.Temperature - block.Target(1);
    if strcmp(string1,'Outlet')
        Out.Measured1 = Inlet.Temperature;
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
end
end%Ends fnction HVAC_NREL