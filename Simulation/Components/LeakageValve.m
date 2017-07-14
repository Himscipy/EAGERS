function Out = LeakageValve(varargin)
%2 inlets: Inlet flow and valve position
%3 outlets: Leakage and Flow1 and flow 2
if length(varargin) == 1%first initialization
    block = varargin{1};
    block.IC = [];
    
    block.InletPorts = {'ValvePos','Intake'};
    block.ValvePos.IC = 0;
    Inlet.ValvePos = block.ValvePos.IC;
    block.Intake.IC.T = 500;
    block.Intake.IC.N2 = 0.0131;
    block.Intake.IC.O2 = 0.0035;
    
    Inlet.Intake = block.Intake.IC;
    spec = fieldnames(Inlet.Intake);
        for i = 1:length(spec)
            Out.Leakage.(spec{i}) = block.leakVal*Inlet.Intake.(spec{i});
            Out.Flow1.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
            Out.Flow2.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*Inlet.ValvePos;
        end
    Out.Leakage.T = Inlet.Intake.T;
    Out.Flow1.T = Inlet.Intake.T;
    Out.Flow2.T = Inlet.Intake.T;

    block.OutletPorts = {'Flow1','Flow2','Leakage'};
    block.Flow2.IC = Out.Flow2;
    block.Leakage.IC = Out.Leakage;
    block.Flow1.IC = Out.Flow1;
    
    block.P_Difference = {};
    Out = block;
elseif length(varargin) == 2
    block = varargin{1};
    Inlet = varargin{2};
    block.ValvePos.IC = Inlet.ValvePos;
    block.Intake.IC = Inlet.Intake;
    spec = fieldnames(Inlet.Intake);
    for i = 1:length(spec)
        Out.Leakage.(spec{i}) = block.leakVal*Inlet.Intake.(spec{i});
        Out.Flow1.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
        Out.Flow2.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*Inlet.ValvePos;
    end
    Out.Leakage.T = Inlet.Intake.T;
    Out.Flow1.T = Inlet.Intake.T;
    Out.Flow2.T = Inlet.Intake.T;
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    if strcmp(string1,'Outlet')
        specName = fieldnames(Inlet.Intake);
        for i = 1:length(specName)
            Out.Leakage.(specName{i}) = block.leakVal*Inlet.Intake.(specName{i});
            Out.Flow1.(specName{i}) = Inlet.Intake.(specName{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
            Out.Flow2.(specName{i}) = Inlet.Intake.(specName{i})*(1-block.leakVal)*Inlet.ValvePos;
        end
        Out.Leakage.T = Inlet.Intake.T;
        Out.Flow1.T = Inlet.Intake.T;
        Out.Flow2.T = Inlet.Intake.T;
    elseif strcmp(string1,'dY')
        %no states
    end
end
end%Ends function LeakageValve