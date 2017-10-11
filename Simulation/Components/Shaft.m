function Out = Shaft(varargin)
% a simple shaft model
% Three (3) inlets: WTurbine, WCompressor, and Gen_Power
% Two (2) outlets: RPM, Steady_Power
% One (1) states: Shaft Speed
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.PMoI = pi()/2*block.Radius^4;%Shaft Polar Moment of Inertia, relates shaft speed to work

    if ~isempty(block.connections{1}) && ~strncmp('Tags',block.connections{1},4)
        r = strfind(block.connections{1},'.');
        block.Scale = ComponentProperty(strcat(block.connections{1}(1:r-1),'.RPMdesign'))/60*(2*pi); %shaft speed in Rad/s normalized by design shaft speed.
    else
        if ~isempty(block.connections{2}) && ~strncmp('Tags',block.connections{2},4)
            r = strfind(block.connections{2},'.');
            block.Scale = ComponentProperty(strcat(block.connections{2}(1:r-1),'.RPMdesign'))/60*(2*pi); %shaft speed in Rad/s normalized by design shaft speed.
        else block.Scale = (block.RPMinit/60*(2*pi));
        end
    end
    block.IC = (block.RPMinit/60*(2*pi))/block.Scale; %shaft speed in Rad/s
    block.UpperBound = 2*block.Scale;
    block.LowerBound = .3*block.Scale;
    
    block.InletPorts = {'WTurbine','WCompressor','Gen_Power'};
    block.WTurbine.IC = 200;%in KW
    block.WTurbine.Saturation = [0,inf];
    block.WCompressor.IC = 100;%in KW
    block.WCompressor.Saturation = [0,inf];
    block.Gen_Power.IC = 100;%in KW
    block.Gen_Power.Saturation = [0,inf];

    block.OutletPorts = {'RPM','Steady_Power'};
    block.RPM.IC = block.RPMinit; %shaft speed in RPM
    block.Steady_Power.IC = 100;%in KW
    
    block.P_Difference = {};
    %no dMdP or mFlow (no pressure calculations)
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    Inlet = checkSaturation(Inlet,block);
    block.Steady_Power.IC = Inlet.WTurbine - Inlet.WCompressor;%in KW
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    if strcmp(string1,'Outlet')
        Out.Steady_Power = Inlet.WTurbine - Inlet.WCompressor;
        Out.RPM = Y(1)*60/(2*pi);%Converts from Radians per Second to RPM
        Tags.(block.name).RPM = Out.RPM;
        Tags.(block.name).SteadyPower = Out.Steady_Power;
    elseif strcmp(string1,'dY')
        dY = 0*Y;
        ShaftPower = (Inlet.WTurbine - Inlet.WCompressor - Inlet.Gen_Power)*1000;%all units should be W
        Moment_Inertia = block.Density * block.Length * block.PMoI;%Moment of Inertia for the shaft
        dY(1) = ShaftPower/(Moment_Inertia*Y(1)); %dw/dt = P/(w*l) units of rad/s
        Out = dY;
    end
end
end%Ends function Shaft