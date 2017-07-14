function Air = makeAir(Tdb,Humidity,massFlow,mode)
%create a flow structure representing air at the correct temperature (C),
%relative humidity(%) and mass flow of dry air (kg/s)
%uses a curve fit for saturation pressure of water in air to calculate
%molar flow of major species (H2O, N2, O2, AR, CO2)
Air.T = Tdb+273.15;
P = 101.325; %kPa, atmospheric pressure
switch mode
    case 'abs' %absolute humidity
        P_H2O = P*Humidity/(0.621945+Humidity);
    case 'rel' %relative humidity
        satP = exp((-5.8002206e3)./Air.T + 1.3914993 - 4.8640239e-2*Air.T + 4.1764768e-5*Air.T.^2 - 1.4452093e-8*Air.T.^3 + 6.5459673*log(Air.T))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
        P_H2O = Humidity/100.*satP;
    case 'dp' %dewpoint
        DP = Humidity;
        satP = exp((-5.8002206e3)./Air.T + 1.3914993 - 4.8640239e-2*Air.T + 4.1764768e-5*Air.T.^2 - 1.4452093e-8*Air.T.^3 + 6.5459673*log(Air.T))/1000; %saturated water vapor pressure ASHRAE 2013 fundamentals eq. 6 in kPa valid for 0 to 200C
        P_H2O = exp(5423*(1/273-1./(273+DP)))./exp(5423*(1/273-1./(273+Air.T))).*satP;%Clausius-Clapeyron equation to calculate relative humidity from dewpoint temperature
end
daMMass = 0.2095*32+0.7809*28+.0092*40+.0004*44;%dry air molar mass
molFlow = massFlow/daMMass./(1-P_H2O/P);%molar flow rate of humid air

Air.H2O = P_H2O/P.*molFlow;
Air.O2 = 0.2095*(1-P_H2O/P).*molFlow;
Air.N2 = 0.7809*(1-P_H2O/P).*molFlow;
Air.AR = 0.0092*(1-P_H2O/P).*molFlow;
Air.CO2 = .0004*(1-P_H2O/P).*molFlow;
end%Ends function makeAir        