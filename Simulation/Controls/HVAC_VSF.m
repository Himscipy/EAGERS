function Out = HVAC_VSF(t,Y, Inlet,block,string1)
% Controls for building HVAC with variable speed fan and proportional cooling
% Three (3) inlets: Building temperature, building humidity, building mode (heating or cooling)
% Three (3) outlets: Supply flow to the building, treated air temperature (T after cooling, T after heating -- different when dehumidifying), Damper position
% One (1) state: mass flow
global Tags
Tol = block.Tolerance;
Terror = Inlet.Temperature - block.Target(1);

if strcmp(string1,'Outlet')
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
        Out.T = {block.Target(2),T};
    else
        Out.T = {T,T};
    end
    Out.Flow = max(min(Y(1)+block.PropGain(1)*Terror,block.maxFlow),block.minFlow);
    Out.Damper = block.Damper.IC ;
    Tags.(block.name).FlowRate = Out.Flow;
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