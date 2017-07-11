function Out = WaterPump(t,Y, Inlet,block,string1)
% Generic water pump for chiller plant
% Two (2) Inlets: Power, Temperature
% One (1) Outlet: water flow
% Zero (0) States: 
global Tags    
if strcmp(string1,'Outlet')
    Flow = Inlet.Power*block.Efficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
    Out.Flow.T = Inlet.Temperature+273.15;
    Out.Flow.H2O = Flow/18;%flow in kmol/s
    Tags.(block.name).Flow_GPM = Out.Flow.H2O*284.96;%flow converted from kmol/s to gallons per minute
elseif strcmp(string1,'dY')  
    %no states
end