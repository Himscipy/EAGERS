function Out = HVAC_VSF(t,Y, Inlet,block,string1)
% Controls for building HVAC with variable speed fan and proportional cooling
% Five (5) inlets: Building temperature, building humidity, building mode (heating or cooling), ambient temperature, ambient humidity
% Three (3) outlets: Supply flow to the building, Cooling load, heating load
% One (1) state: mass flow
global Tags
Tol = block.Tolerance;
Terror = Inlet.Temperature - block.Target(1);
flow = max(min(Y(1)+block.PropGain(1)*Terror,block.maxFlow),block.minFlow);
%ambient air
AmbientAir = makeAir(Inlet.Tamb,Inlet.ambHumidity,(1-block.Entrainment/100)*flow,'abs');
%recirculated air
RecircAir = makeAir(Inlet.Temperature,Inlet.Humidity,block.Entrainment/100*flow,'abs');
%mixed air
MixedAir = MixAir(RecircAir,AmbientAir);
mixed_AH = MixedAir.H2O*18/(MassFlow(MixedAir)-MixedAir.H2O*18);%actual humidity: kg H20/kg dry air

if strcmp(string1,'Outlet')
    if strcmp(Inlet.Mode,'cooling')
        if Terror<-Tol ||  (Tags.(block.name).Cooling<1e-6 && Terror<Tol) %don't run AC
            T = MixedAir.T;
        elseif mixed_AH>block.maxHumidity && Terror<0%if de-humidifying, and too cold (can't shut off cooling, so re-heat extra)
            T = (block.ColdAirSetpoint + block.Target(1))/2 + 273.15;
        else %if Terror>2 || (Tags.(block.name).Cooling>0 && Terror>-2)
            T = block.ColdAirSetpoint+273.15;
        end
    else
        if Terror>Tol ||  (Tags.(block.name).Heating<1e-6 && Terror>-Tol) %don't run heater
            T = MixedAir.T;
        else %if Terror<-2 || (Tags.(block.name).Heating>0 && Terror<2)
            T = block.HotAirSetpoint+273.15;
        end
    end
    
    %%cooling
    if mixed_AH>block.maxHumidity %if de-humidifying, cool to the dew point
        CooledAir = makeAir(block.Target(2),100,flow,'rel');
    elseif T<MixedAir.T % cool to the supply temperature if cooling (no re-heating
        CooledAir = MixedAir;
        CooledAir.T = T;%cooled mixed air
    else %no cooling (only heating)
        CooledAir = MixedAir;
    end
    Tags.(block.name).Cooling  = enthalpyAir(MixedAir)*(MassFlow(MixedAir) - MixedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    %%heating
    HeatedAir = CooledAir;
    HeatedAir.T = T;
    Tags.(block.name).Heating  = enthalpyAir(HeatedAir)*(MassFlow(HeatedAir) - HeatedAir.H2O*18) - enthalpyAir(CooledAir)*(MassFlow(CooledAir) - CooledAir.H2O*18);
    
    Out.Supply = HeatedAir;
    Tags.(block.name).FlowRate = MassFlow(Out.Supply) - Out.Supply.H2O*18;
    Tags.(block.name).Temperature = Out.Supply.T - 273.15;
    Tags.(block.name).Humidity = Out.Supply.H2O*18/Tags.(block.name).FlowRate; % kg H2O / kg dry air
    Tags.(block.name).CoolingPower = Tags.(block.name).Cooling/block.COP;
elseif strcmp(string1,'dY')
    dY = Y*0;
    %%mass flow
    if Tags.(block.name).Cooling == 0 % if not cooling, run at minimum flow
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