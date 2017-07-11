function block = InitializeHVAC_VSF(varargin)
% Controls for building HVAC with variable speed fan and proportional cooling
% Three (3) inlets: Building temperature, building humidity, building mode (heating or cooling)
% Three (3) outlets: Supply flow to the building, treated air temperature (T after cooling, T after heating -- different when dehumidifying), Damper position
% One (1) state: mass flow
global Tags
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Mass flow (of dry air)';};
    s = block.connections{1};
    r = strfind(s,'.');
    block.building = s(1:r(1)-1);
    
    block.InletPorts = {'Temperature','Humidity','Mode'};
    block.Temperature.IC = 22.2; 
    block.Humidity.IC = 0.0085; 
    block.Mode.IC = 'cooling';
        
    block.OutletPorts = {'Flow';'Tset';'Damper'};   
    block.Flow.IC = block.minFlow;
    block.Tset.IC = {block.ColdAirSetpoint;block.ColdAirSetpoint;};
    block.Damper.IC = (1-block.Entrainment/100);
    
    block.InitializeError = 1;
    
    block.Scale = [block.maxFlow;];
    block.IC = [block.minFlow/block.maxFlow;]; % inital condition 
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    %during initialization we are providing the cooling or heating which
    %results in steady-state temperature at the target set-point
    %the controller may cycle the cooling
    Terror = Inlet.Temperature - block.Target(1);
    %find temperature and flow
    if strcmp(Inlet.Mode,'cooling')
        flow = block.Flow.IC; %mass flow of dry air
        if Terror<0 && flow-block.minFlow<1e-10 %if its too cold and the flow has already been reduced, raise the supply temperature
            T = block.Tset.IC{2} - Terror;
        elseif Terror>0 && block.maxFlow-flow<1e-10 %too hot & max flow
            T = block.Tset.IC{2} - Terror;
        else
            T = block.ColdAirSetpoint;
            block.Flow.IC = min(block.maxFlow,max(block.minFlow,(1 + Terror/5)*flow));
        end
    else %heating
        block.Flow.IC = block.minFlow;
        T = block.Tset.IC{2} - Terror;
    end
    if Inlet.Humidity>block.maxHumidity %if de-humidifying, cool to the dew point
        block.Tset.IC = {block.Target(2),T};
    else
        block.Tset.IC = {T,T};
    end
    block.IC = [flow/block.maxFlow;];
    block.InitializeError = abs(Terror);
    Tags.(block.name).CoolingPower = Tags.(block.building).Cooling/block.COP;
end