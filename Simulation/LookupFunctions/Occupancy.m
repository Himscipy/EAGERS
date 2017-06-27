function Out = Occupancy(t)
%Function outputting the occupancy given the schedule for a Large Office 

AreaLoad = 5.38195521e-2;%W/m^2
A = 313.42;%m^2
MaxLoad = AreaLoad*A;
Fraction =[0.05 0 0 0 0 0 0 .1 .2 .95 .95 .95 .95 .5 .95 .95 .95 .95 .3 .1 .1 .1 .1 .05 .05];
Load = Fraction*MaxLoad;
hour = linspace(0,24,25);
h = mod(t/3600,24);
Out = interp1(hour,Load,h);
