function Out = SimpleMix(varargin) 
% a simple adiabatic mixing model mixing two or more inlet flows 
% Multiple (n+1) inlets: Flows (consisting of Temperature and flow rates of individual species) and P out 
% Three (3) outlets: Pressure, Flow , and Temperature 
% Two (2) states: Temperature and pressure
global Tags
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.Pdrop = 2; % presure drop across mixing volume (kPa)
    block.Scale = zeros(2,1);
    block.IC = ones(length(block.Scale),1);
    block.UpperBound = [inf,inf];
    block.LowerBound = [0,0];
    
    block.Scale(1) = block.Tinit; 
    block.Scale(2) = 101 + block.Pdrop;
   
    block.InletPorts = {};
    for i = 1:1:block.inlets
        name = strcat('Inlet',num2str(i));
        block.(name) = [];
        block.InletPorts(end+1) = cellstr(name);
        block.(name).Saturation = [0,inf];
    end
    
    block.InletPorts(end+1) = {'Pout'};
    block.Pout.IC = 101;
    block.Pout.Saturation = [0,inf];
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'Outlet';'Temperature';'Pin';};
    block.Outlet.IC = block.SpeciesInit;
    block.Outlet.IC.T = block.Tinit;
    block.Temperature.IC = block.Outlet.IC.T;
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.P_Difference = {'Pin','Pout'};
    %no dMdP or mFlow (fixed pressure drop)
    Out = block;
elseif length(varargin)==2%% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet =varargin{2};
    Inlet = checkSaturation(Inlet,block);
    n = block.inlets;     

    inlets = fieldnames(Inlet);
    spec = {};
    for i = 1:1:n
        spec2 = fieldnames(Inlet.(inlets{i}));
        spec2 = spec2(~strcmp('T',spec2));
        for j = 1:1:length(spec2)
            if ismember(spec2{j},spec)
                NetIn.(spec2{j}) = NetIn.(spec2{j}) + Inlet.(inlets{i}).(spec2{j});
            else
                spec(end+1) = spec2(j);
                NetIn.(spec2{j}) = Inlet.(inlets{i}).(spec2{j});
            end
        end
    end
    
    NetFlowIn = NetFlow(NetIn);
    block.Pfactor = NetFlowIn/block.Pdrop;
    
    H = zeros(n,1);
    Tin = zeros(n,1);
    for i = 1:1:n
        [~,H(i)] = enthalpy(Inlet.(inlets{i}));
        Tin(i) = Inlet.(inlets{i}).T;
    end
    
    Hin = sum(H);
    NetOut = NetIn;
    NetOut.T = mean(Tin);
    Terror = 1;
    while abs(Terror)>.01
        [~,Hout] = enthalpy(NetOut);
        Cp = SpecHeat(NetOut);
        Terror = (Hin-Hout)/(Cp*NetFlowIn);
        NetOut.T = NetOut.T+Terror;
    end
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.Outlet.IC = NetOut; 
    block.Temperature.IC = NetOut.T; 
    
    block.Scale = zeros(2,1);
    block.IC = [1,1];
    block.Scale = [NetOut.T;block.Pin.IC;];
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Inlet = checkSaturation(Inlet,block);
    Pin =Y(end);
    NetOut.T = Y(1);
    inlets = fieldnames(Inlet);
    n = block.inlets;
    H = zeros(n,1);
    spec ={};
    for i = 1:1:n
        spec2 = fieldnames(Inlet.(inlets{i}));
        spec2 = spec2(~strcmp('T',spec2));
        for j = 1:1:length(spec2)
            if ~ismember(spec2{j},spec)
                spec(end+1) = spec2(j);
                NetIn.(spec2{j}) = Inlet.(inlets{i}).(spec2{j});
            else
                NetIn.(spec2{j}) = NetIn.(spec2{j}) + Inlet.(inlets{i}).(spec2{j});
            end
        end
        H(i) = enthalpy(Inlet.(inlets{i}));
    end

    scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/NetFlow(NetIn);%total cold flow out
    for i = 1:1:length(spec)
        NetOut.(spec{i}) = NetIn.(spec{i})*scaleFlow;
    end

    if strcmp(string1,'Outlet')
        Out.Pin = Pin;
        Out.Outlet = NetOut;
        Out.Temperature = Out.Outlet.T;
        Tags.(block.name).MassFlow = MassFlow(NetOut);
        Tags.(block.name).Temperature = NetOut.T;
    elseif strcmp(string1,'dY')
        dY = 0*Y;
        Hin = sum(H);  
        %temperature
        scaleFlow2 = NetFlow(NetIn)/NetFlow(NetOut);
        Hout = enthalpy(NetOut)*scaleFlow2; %scale Hin and Hout to the same flow rate (no reactions so easy)
        Cp = SpecHeat(NetOut);
        dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*NetOut.T*block.Ru;
        %pressure
        dY(end) = (NetFlow(NetIn) - NetFlow(NetOut))*block.Ru*Y(1)/(block.Vol);
        Out = dY;
    end
end
end%Ends function SimpleMix
