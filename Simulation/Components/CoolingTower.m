function Out = CoolingTower(t,Y, Inlet,block,string1)
% Water cooling tower
% Four (4) Inlets: supply water flow, Fan power, Ambient temperature, ambient humididty
% Teo (2) Outlets: Chilled return water flow, Chilled water temperature
% Three (3) States: Air temperature after cooling, heat exchanger temperature, chilled water temperature
global Tags    
if strcmp(string1,'Outlet')
    Out.Return.H2O = Inlet.Supply.H2O;
    Out.Return.T = Y(3)+273.15;
    Out.Temperature = Y(3);%temperature in Celcius
    Tags.(block.name).Temperature = Y(3);
elseif strcmp(string1,'dY')
    Cp_H2O = 4.186;
    Cp_Air = 1.006;
    P = 101.325;
    P_H2O = P*Inlet.ambHumidity/(0.621945+Inlet.ambHumidity);
    Tamb_K = Inlet.Tamb + 273.15;
    satP = exp((-5.8002206e3)./Tamb_K + 1.3914993 - 4.8640239e-2*Tamb_K + 4.1764768e-5*Tamb_K.^2 - 1.4452093e-8*Tamb_K.^3 + 6.5459673*log(Tamb_K))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
    RH = P_H2O/satP*100;
    Tair_WB = Inlet.Tamb*atan(0.151977*(RH + 8.313659)^.5) + atan(Inlet.Tamb + RH) - atan(RH - 1.676331) + 0.00391838*(RH)^(3/2)*atan(0.023101*RH) - 4.686035; % http://journals.ametsoc.org/doi/full/10.1175/JAMC-D-11-0143.1
    AirFlow = Inlet.Power*block.FanEfficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
    %co-flow
    Qhot = block.hconv*block.Area*(Y(2)-Y(1)); %heat transfer from air to heat exchanger
    Qcold = block.hconv*block.Area*(Y(3)-Y(2)); %heat transfer from heat exchanger to water
%     %counterflow
%     Qhot = block.hconv*block.Area*((Y(1)+Tair_WB)/2-Y(2)); %heat transfer from air to heat exchanger
%     Qcold = block.hconv*block.Area*(Y(2) - (Y(3)+(Inlet.Supply.T-273.15))/2); %heat transfer from heat exchanger to water
    
    dY(1) = ((Tair_WB + Qhot/(Cp_Air*AirFlow))-Y(1))/1;
    dY(2) = (Qcold - Qhot)/(block.CP*block.Mass);
    dY(3) = ((Inlet.Supply.T-273.15) - Qcold/(Cp_H2O*Inlet.Supply.H2O*18) - Y(3))/10;
    Out = dY;
    Tags.(block.name).WetBulbTemperature = Tair_WB;
end