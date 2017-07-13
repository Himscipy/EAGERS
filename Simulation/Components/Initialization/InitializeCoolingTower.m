function block = InitializeCoolingTower(varargin)
% Generic Cooling Tower for chiller plant
% Four (4) Inlets: Cold water flow, Fan power, Ambient temperature, ambient humidity
% Two (2) Outlets: Chilled water fow, Chilled water temp
% Three (3) States: WB air temp out, HX temp, chilled water temp (all in Celcius)
global Tags
block = varargin{1};
Cp_H2O = 4.186; %kJ/kg
if length(varargin)==1 % first initialization
    block.Area = 2*block.Cooling_Tons*3.517/(block.hconv*block.HeaterDeltaT); %Surface area (m^2)
    
    block.Scale = [20 21 22]; %temperatures
    block.IC = ones(length(block.Scale),1);
    %%
    block.InletPorts = {'Power','Supply','Tamb','ambHumidity'};
    block.Power.IC = 1;
    block.Supply.IC.T = 25+273.15;
    block.Supply.IC.H2O = block.Cooling_Tons*3.517/(Cp_H2O*(block.Supply.IC.T-273.15-block.Scale(3)));
    block.Tamb.IC = 25;
    block.ambHumidity.IC = 0.01;
    
    block.OutletPorts = {'Return';'Temperature';};
    block.Return.IC.T = 22+273.15;
    block.Return.IC.H2O = block.Supply.IC.H2O;
    block.Temperature.IC = block.Return.IC.T-273.15;
    
    block.P_Difference = {};
end
if length(varargin)==2 %% Have inlets connected, re-initialize  
    Inlet = varargin{2};
    Cp_Air = 1.006; %kJ/kmol
%     rho_air = 1.204; %density of dry air at sea level
    AirFlow = Inlet.Power*block.FanEfficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
    P = 101.325;
    P_H2O = P*Inlet.ambHumidity/(0.621945+Inlet.ambHumidity);
    Tamb_K = (Inlet.Tamb+273.15);
    satP = exp((-5.8002206e3)./Tamb_K + 1.3914993 - 4.8640239e-2*Tamb_K + 4.1764768e-5*Tamb_K.^2 - 1.4452093e-8*Tamb_K.^3 + 6.5459673*log(Tamb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
    RH = P_H2O/satP*100;
    Tair_WB = Inlet.Tamb*atan(0.151977*(RH + 8.313659)^.5) + atan(Inlet.Tamb + RH) - atan(RH - 1.676331) + 0.00391838*(RH)^(3/2)*atan(0.023101*RH) - 4.686035; % http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
    error = 1;
    
    a = .75; %change at most a degrees C each iteration
    while sum(abs(error))>1e-3
        dT = (block.Scale(3) - block.Scale(1));
        Q_HT = block.hconv*block.Area/2*dT; %coflow
        Q_C = Cp_H2O*Inlet.Supply.H2O*18*(Inlet.Supply.T-273.15 - block.Scale(3)); %how much cooling to reach previous guess of return water temp
        dT = dT*(1+.5*(Q_C/Q_HT-1));
        Q = (Q_HT + Q_C)/2;
        
%         Q = block.hconv*block.Area/2*((block.Scale(3)+Inlet.Supply.T-273.15)/2 - (Tair_WB+block.Scale(1))/2); %counterflow
        error(1) = (Tair_WB + Q/(Cp_Air*AirFlow)) - block.Scale(1);
        error(2) = (Inlet.Supply.T-273.15 - Q/(Cp_H2O*Inlet.Supply.H2O*18)) - block.Scale(3);
        block.Scale(1) = block.Scale(1) + max(-a,min(a,error(1)));
        block.Scale(3) = block.Scale(1) + dT;
    end
    block.Scale(2) = (block.Scale(1) + block.Scale(3))/2;
    
    block.Return.IC.T = block.Scale(3) + 273.15;
    block.Return.IC.H2O = Inlet.Supply.H2O;
    block.Temperature.IC = block.Scale(3);
    Tags.(block.name).Temperature = block.Temperature.IC;
    Tags.(block.name).WetBulbTemperature = Tair_WB;
end