function Out = InternalHeatGains(t)
%Internal Gains based off e+ Large Office
A = 313.42;%m^2
hour = linspace(1,24,24);

%% Lighting Heat Gain
LightSched = [.04619 .05 .05 .05 .05 .05 .1 .09238 .27714 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .4619 .27714 .27714 .18476 .18476 .09238 .04619];
MaxLightLoad = 9.6875*A;
LightLoad = LightSched*MaxLightLoad;

%% Solar Heat Gain
SolarLoad = 0;%Unknown Value

%% Plug Load Heat Gain
%Loads will be different depending on zone.  Data Centers are always on with higher plug loads
PlugSched = [.4 .4 .4 .4 .4 .4 .4 .4 .4 .9 .9 .9 .9 .8 .9 .9 .9 .9 .5 .4 .4 .4 .4 .4 .4];
MaxPlug = 8.07007*A;
PlugLoad = PlugSched*MaxPlug;

%% Total Internal Gain

h = mod(t/3600,24);
Out = interp1(hour,LightLoad,h) + SolarLoad + interp1(hour,PlugLoad,h);