function block = InitializeWaterPump(varargin)
% Generic water pump for chiller plant
% Two (2) Inlets: Power, Temperature
% One (1) Outlet: water flow
% Zero (0) States: 
global Tags
block = varargin{1};
if length(varargin)==1 % first initialization
    block.Scale = []; 
    block.IC = [];
    %%
    block.InletPorts = {'Power','Temperature'};
    block.Power.IC = 1;
    block.Temperature.IC = 12+273.15;
    
    block.OutletPorts = {'Flow';};
    block.Flow.IC.T = block.Temperature.IC ;
    block.Flow.IC.H2O = 1;
    
    block.P_Difference = {};
    
    Tags.(block.name).Flow_GPM = block.Flow.IC.H2O*284.96;%flow converted from kmol/s to gallons per minute
end
if length(varargin)==2 %% Have inlets connected, re-initialize    
    Inlet = varargin{2};
    Flow = Inlet.Power*block.Efficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
    block.Flow.IC.T = Inlet.Temperature+273.15;
    block.Flow.IC.H2O = Flow/18;%flow in kmol/s
    Tags.(block.name).Flow_GPS = block.Flow.IC.H2O*284.96;%flow converted from kmol/s to gallons per minute
end