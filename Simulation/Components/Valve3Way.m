function Out = Valve3Way(varargin)
% a simple flow splitting block with 1 inlet flow,and 2 outlets
% Two (2) inlets: inlet flow ,  valve postion
% Two (2) outlets:  Flow1 , and Flow2 
% Zero (0) states:
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.IC =[];% no states

    if ischar(block.InitialFlowIn)
        block.InitialFlowIn = ComponentProperty(block.InitialFlowIn);
    end
    spec = fieldnames(block.InitialFlowIn);
    for i = 1:1:length(spec)
        FlowIn.(spec{i}) = block.InitialFlowIn.(spec{i});
        Out1.(spec{i}) = block.InitialFlowIn.(spec{i})*block.PercOpen;
        Out2.(spec{i}) = block.InitialFlowIn.(spec{i})*(1-block.PercOpen);
    end
    FlowIn.T = block.InitialFlowIn.T;
    Out1.T = block.InitialFlowIn.T;
    Out2.T = block.InitialFlowIn.T;
    
    block.InletPorts = {'Inlet','ValvePos'};
    block.Inlet.IC = FlowIn;
    block.Inlet.Saturation = [0,inf];
    block.ValvePos.IC = block.PercOpen;
    block.ValvePos.Saturation = [0,1];
    
    block.OutletPorts = {'Out1','Out2'};
    block.Out1.IC = Out1;
    block.Out2.IC = Out2;  
    
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    spec = fieldnames(Inlet.Inlet);
    for i = 1:1:length(spec)
        if ~strcmp(spec{i},'T')
            Out1.(spec{i}) = Inlet.Inlet.(spec{i})*Inlet.ValvePos;
            Out2.(spec{i}) = Inlet.Inlet.(spec{i})*(1-Inlet.ValvePos);
        end
    end
    Out1.T = Inlet.Inlet.T;
    Out2.T = Inlet.Inlet.T;
    block.Out1.IC = Out1;
    block.Out2.IC = Out2;  
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    if strcmp(string1,'Outlet')
        spec = fieldnames(Inlet.Inlet);
        for i = 1:length(spec)
            Out.Out1.(spec{i}) = Inlet.ValvePos*Inlet.Inlet.(spec{i});
            Out.Out2.(spec{i}) = (1-Inlet.ValvePos)*Inlet.Inlet.(spec{i});
        end
        Out.Out1.T = Inlet.Inlet.T;
        Out.Out2.T = Inlet.Inlet.T;
    elseif strcmp(string1,'dY')
        %no states
    end
end
end%Ends function Valve3Way