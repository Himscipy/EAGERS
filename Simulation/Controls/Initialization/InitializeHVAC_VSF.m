function block = InitializeHVAC_VSF(varargin)
% Controls for building HVAC with variable speed fan and proportional cooling
% Five (5) inlets: Building temperature, building humidity, building mode (heating or cooling), ambient temperature, ambient humidity
% Three (3) outlets: Supply flow to the building, Cooling load, heating load
% One (1) state: mass flow
global Tags
block = varargin{1};
P = 101.325; %kPa, atmospheric pressure value
if length(varargin)==1 %first initialization
    block.description = {'Air Supply to Building';};

    block.InletPorts = {'Temperature','Humidity','Mode','Tamb','ambHumidity'};
    block.Temperature.IC = 22.2; 
    block.Humidity.IC = 0.0085; 
    block.Mode.IC = 'cooling';
    block.Tamb.IC = 25;
    block.ambHumidity.IC = .01;
        
    block.OutletPorts = {'Supply';'Cooling';'Heating';};   
    block.Supply.IC = makeAir(block.ColdAirSetpoint,50,block.minFlow,'rel');
    block.Cooling.IC = 0;
    block.Heating.IC = 0;
    
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
        flow = MassFlow(block.Supply.IC) - block.Supply.IC.H2O*18; %mass flow of dry air
        if Terror<0 && flow-block.minFlow<1e-10 %if its too cold and the flow has already been reduced, raise the supply temperature
            T = block.Supply.IC.T - Terror;
        elseif Terror>0 && block.maxFlow-flow<1e-10 %too hot & max flow
            T = block.Supply.IC.T - Terror;
        else
            T = block.ColdAirSetpoint + 273.15;
            flow = min(block.maxFlow,max(block.minFlow,(1 + Terror/5)*flow));
        end
    else %heating
        flow = block.minFlow;
        T = block.Supply.IC.T - Terror;
    end

    %Mix fresh and recirculated air at this flow rate
    AmbientAir = makeAir(Inlet.Tamb,Inlet.ambHumidity,(1-block.Entrainment/100)*flow,'abs');
    %recirculated air
    RecircAir = makeAir(Inlet.Temperature,Inlet.Humidity,block.Entrainment/100*flow,'abs');
    %mixed air
    MixedAir = MixAir(RecircAir,AmbientAir);
    mixed_AH = MixedAir.H2O*18/(MassFlow(MixedAir)-MixedAir.H2O*18);%actual humidity: kg H20/kg dry air
    
    %%cooling
    if mixed_AH>block.maxHumidity %if de-humidifying, cool to the dew point
        CooledAir = makeAir(block.Target(2),100,flow,'rel');%saturated air
        block.Cooling.IC = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    elseif (T-273.15)<block.Target(1) % cool to the supply temperature if cooling
        CooledAir = MixedAir;
        CooledAir.T = T;%cooled mixed air
        block.Cooling.IC = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    else
        CooledAir = MixedAir;
        block.Cooling.IC = 0;
    end
    %heating
    HeatedAir = CooledAir;
    HeatedAir.T = T;
    block.Heating.IC = enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    
    block.Supply.IC = HeatedAir;
    
    block.IC = [flow/block.maxFlow;];
    block.InitializeError = abs(Terror); %0.5*block.InitializeError + 0.5*
    Tags.(block.name).Cooling = block.Cooling.IC;
    Tags.(block.name).Heating = block.Heating.IC;
    Tags.(block.name).CoolingPower = block.Cooling.IC/block.COP;
end