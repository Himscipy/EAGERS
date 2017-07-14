function Out = Source(varargin)
% a source of regular fuel with no states
% the type of fuel (species) must be initialized first
% Two (2) inlets: Temperature and flow rate
% no state (will only ever get called upon to calculate outlet)
global Tags
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.IC = [];%no states exist
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Temperature','Species','Flow'};
    block.Temperature.IC = 298; 
    block.Species.IC = block.InitialComposition; 
    block.Flow.IC = 1;  
    
    block.OutletPorts = {'Outlet'};
    block.Outlet.IC.T = block.Temperature.IC;
    speciesNames = fieldnames(block.Species.IC);
    for i = 1:1:length(speciesNames)
        block.Outlet.IC.(speciesNames{i}) = block.Species.IC.(speciesNames{i})*block.Flow.IC;
    end
    
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    block.Temperature.IC = Inlet.Temperature;
    block.Species.IC = Inlet.Species;  
    block.Flow.IC = Inlet.Flow;
    block.Outlet.IC.T = Inlet.Temperature;
    speciesNames = fieldnames(Inlet.Species);
    for i = 1:1:length(speciesNames)
        block.Outlet.IC.(speciesNames{i}) = Inlet.Species.(speciesNames{i})*Inlet.Flow;
    end
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    if strcmp(string1,'Outlet')
        Out.Outlet.T = Inlet.Temperature;
        speciesNames = fieldnames(Inlet.Species);
        for i = 1:1:length(speciesNames)
            Out.Outlet.(speciesNames{i}) = Inlet.Species.(speciesNames{i})*Inlet.Flow;
        end
        Tags.(block.name).Outlet = Out.Outlet;
    elseif strcmp(string1,'dY')
        %no states
    end
end
end%Ends function Source