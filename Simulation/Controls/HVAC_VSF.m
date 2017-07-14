function Out = HVAC_VSF(varargin)
% Controls for building HVAC with variable speed fan and proportional cooling
% Three (3) inlets: Building temperature, building humidity, building mode (heating or cooling)
% Three (3) outlets: Supply flow to the building, treated air temperature (T after cooling, T after heating -- different when dehumidifying), Damper position
% One (1) state: mass flow
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
    block.InletPorts(end+1:end+3) = {'Temperature','Humidity','Mode'};
    block.Temperature.IC = 22.2;
    block.Humidity.IC = 0.0085; 
    block.Mode.IC = 'cooling';
    
    block.OutletPorts(end+1:end+3) = {'Flow';'Tset';'Damper'}; 
    block.Flow.IC = block.minFlow;
    block.Tset.IC = {block.ColdAirSetpoint;block.ColdAirSetpoint;};
    block.Damper.IC = (1-block.Entrainment/100);  
    
    s = block.connections{3};
    r = strfind(s,'.');
    block.building = s(1:r(1)-1);
    
    block.Scale = [block.maxFlow;];
    block.IC = [block.minFlow/block.maxFlow;]; % inital condition 
    block.UpperBound = inf;
    block.LowerBound = 0;
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
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
    block.Measured1.IC = Inlet.Temperature;
    P = 101.325;
    P_H2O = P*Inlet.Humidity/(0.621945+Inlet.Humidity);
    T_K = Inlet.Temperature;
    satP = exp((-5.8002206e3)./T_K + 1.3914993 - 4.8640239e-2*T_K + 4.1764768e-5*T_K.^2 - 1.4452093e-8*T_K.^3 + 6.5459673*log(T_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
    RH = P_H2O/satP*100;
    T_WB = Inlet.Temperature*atan(0.151977*(RH + 8.313659)^.5) + atan(Inlet.Temperature + RH) - atan(RH - 1.676331) + 0.00391838*(RH)^(3/2)*atan(0.023101*RH) - 4.686035; % http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
    block.Measured2.IC = T_WB;
    
    block.IC = [flow/block.maxFlow;];
    block.InitializeError = abs(Terror);
    Tags.(block.name).CoolingPower = Tags.(block.building).Cooling/block.COP;
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Tol = block.Tolerance;
    Terror = Inlet.Temperature - block.Target(1);

    if strcmp(string1,'Outlet')
        Out.Measured1 = Inlet.Temperature;
        if strcmp(Inlet.Mode,'cooling')
            if Terror<-Tol ||  (Tags.(block.building).Cooling<1e-6 && Terror<Tol) %don't run AC
                T = [];
            else %if Terror>2 || (Tags.(block.name).Cooling>0 && Terror>-2)
                T = block.ColdAirSetpoint;
            end
        else
            if Terror>Tol ||  (Tags.(block.building).Heating<1e-6 && Terror>-Tol) %don't run heater
                T = [];
            else %if Terror<-2 || (Tags.(block.name).Heating>0 && Terror<2)
                T = block.HotAirSetpoint;
            end
        end
        if Inlet.Humidity>block.maxHumidity %if de-humidifying, cool to the dew point
            Out.Tset = {block.Target(2),T};
        else
            Out.Tset = {T,T};
        end
        Out.Flow = max(min(Y(1)+block.PropGain(1)*Terror,block.maxFlow),block.minFlow);
        Out.Damper = block.Damper.IC ;
        Tags.(block.name).FlowRate = Out.Flow;
        if isempty(T)
            T = Inlet.Temperature;
        end
        Tags.(block.name).Temperature = T;
        Tags.(block.name).CoolingPower = Tags.(block.building).Cooling/block.COP;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        %%mass flow
        if Tags.(block.building).Cooling == 0 % if not cooling, run at minimum flow
            dY(1) = (block.minFlow - Y(1))/60;        
        else%control mass flow to control room cooling
            if Terror<0 && Y(1)<= block.minFlow %room is too cold, and already at minimum flow, maintain same flow rate (eventually cooling will shut off)
                dY(1) = (block.minFlow - Y(1))/600;
            elseif Terror>0 && Y(1)>=block.maxFlow %room is too hot, and already at max flow, maintain same flow rate (maximum cooling reached)
                dY(1) = (block.maxFlow - Y(1))/600;
            else %room is too cold, reduce flow rate, room is too hot, increase flow
                dY(1) = block.Gain(1)*Terror;
            end
        end
        Out = dY;
    end
end
end%Ends function HVAC_VSF