function Out = WaterPump(varargin)
% Generic water pump for chiller plant
% Two (2) Inlets: Power, Temperature
% One (1) Outlet: water flow
% Zero (0) States: 
global Tags    
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.IC = [];
    %%
    block.InletPorts = {'Power','Temperature'};
    block.Power.IC = 1;
    block.Temperature.IC = 12+273.15;
    
    block.OutletPorts = {'Flow';};
    block.Flow.IC.T = block.Temperature.IC ;
    block.Flow.IC.H2O = 1;
    
    block.P_Difference = {};
    
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize    
    block = varargin{1};
    Inlet = varargin{2};
    Flow = Inlet.Power*block.Efficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
    block.Flow.IC.T = Inlet.Temperature+273.15;
    block.Flow.IC.H2O = Flow/18;%flow in kmol/s
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    if strcmp(string1,'Outlet')
        Flow = Inlet.Power*block.Efficiency*1e3/(9.81*block.HeadLoss); % flow rate in kg/s
        Out.Flow.T = Inlet.Temperature+273.15;
        Out.Flow.H2O = Flow/18;%flow in kmol/s
        Tags.(block.name).Flow_GPM = Out.Flow.H2O*284.96;%flow converted from kmol/s to gallons per minute
    elseif strcmp(string1,'dY')  
        %no states
    end
end
end%Ends function WaterPump