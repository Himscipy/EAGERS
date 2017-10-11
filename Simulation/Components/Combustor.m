function Out = Combustor(varargin)
% a combustor with 2 inlet flows (air and fuel) modeled with 3 reactions
% The completeness of each reaction is adjusted with Rxn1 --3;
% The combustor is divided into a core flow and a bypass flow
% Four (4) inlets: air flow (TXN), fuel flow (TXN), bypass split between core and bypass flow, and outlet pressure
% Three (3) outlets: Bypass flow, core flow and pressure at the inlet
% Six (6) states: Tcombustor, Tbypass, Twall, Tcasing, Equivelance ratio Inlet Pressure
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.Scale = [1800; 830; 870; 450; 1; 101;];% Tcombustor, Tbypass, Twall, Tcasing, Equivalance ratio, Pcomb
    block.IC = [ones(4,1); 0.5;1;];
    block.UpperBound = [inf,inf,inf,inf,inf,inf];
    block.LowerBound = [0,0,0,0,0,0];
    
    block.Diameter = 0.2;
    block.Length = 2;
    
    block.Vol = pi/4*block.Diameter^2*block.Length;
    
    block.SpecHeatWall = .5; % Specific Heat of Combustor (kJ/kg*K)
    block.SpecHeatCasing = .5; % Specific Heat of Combustor (kJ/kg*K)
    block.Tamb = 330;
    block.Pdrop = 2; % presure drop across combustor (kPa)
    block.Rxn1 = 1; % CH4+1.5O2 --> CO+ 2H2O
    block.Rxn2 = .7; % CO + .5O2 --> CO2
    block.Rxn3 = 1; % H2 + .5O2 --> H2O
    
    block.SurfA = pi*block.Diameter*block.Length; %combustor surface area
    block.AmbConv = 5; % Ambient Convection Coefficient
    block.AmbEpsilon = .8; %Ambient Radiation Epsilon
    block.AmbSigma = 5.67e-8; %Ambient Radiation Sigma

    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Air','Fuel','Pout'};
    
    block.Air.IC.T = 700;
    block.Air.IC.N2 = .79;
    block.Air.IC.O2 = .21;
    block.Air.Saturation = [0,inf];
    
    block.Fuel.IC.T = 700;
    block.Fuel.IC.CH4 = 0.9;
    block.Fuel.IC.CO = 0.05;
    block.Fuel.IC.CO2 = 0.03;
    block.Fuel.IC.N2 = 0.02;
    block.Fuel.Saturation = [0,inf];
    
    block.Pout.IC = 100;
    block.Pout.Saturation = [0,inf];
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'Bypass','Main','Pin'};
    
    block.Bypass.IC.T = 1000;
    block.Bypass.IC.N2 = 0.79;
    block.Bypass.IC.O2 = 0.21;
    
    block.Main.IC.T = 2000;
    block.Main.IC.H2O = 0.1;
    block.Main.IC.N2 = .75;
    block.Main.IC.O2 = .15;
    
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = 5; %identifies the state # of the pressure state if this block has one

    block.P_Difference = {'Pin','Pout'};
    %no dMdP or mFlow (fixed pressure drop)
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    
    Inlet = checkSaturation(Inlet,block);
    Bypass = 1 - (2*Inlet.Fuel.CH4+.5*Inlet.Fuel.CO+.5*Inlet.Fuel.H2)/(Inlet.Air.O2*block.EquivSet);
    R(1) = Inlet.Fuel.CH4*block.Rxn1;
    R(2) = (Inlet.Fuel.CO+R(1))*block.Rxn2;
    R(3) = Inlet.Fuel.H2*block.Rxn3;
    
    sumR = 1.5*R(1) + 0.5*R(2) + 0.5*R(3);
    
    if sumR>Inlet.Air.O2*(1-Bypass) %rich combustion limited by O2
        partialComb = Inlet.Air.O2*(1-Bypass)/sumR;
        R = partialComb*R;
    end

    specinterest = {'CH4','CO','CO2','H2','H2O','N2','O2'};
    for i =1:1:length(specinterest)
        if ~isfield(Inlet.Air,specinterest{i})
            Inlet.Air.(specinterest{i}) = 0;
        end
        if ~isfield(Inlet.Fuel,specinterest{i})
            Inlet.Fuel.(specinterest{i}) = 0;
        end
    end

    block.Main.IC.CH4 = Inlet.Fuel.CH4-R(1);
    block.Main.IC.CO = Inlet.Fuel.CO + Inlet.Air.CO*(1-Bypass) + R(1)- R(2);
    block.Main.IC.CO2 = Inlet.Fuel.CO2 + Inlet.Air.CO2*(1-Bypass) + R(2);
    block.Main.IC.H2 = Inlet.Fuel.H2 - R(3);
    block.Main.IC.H2O = Inlet.Fuel.H2O + Inlet.Air.H2O*(1-Bypass) + 2*R(1) + R(3);
    block.Main.IC.N2 = Inlet.Fuel.N2 + Inlet.Air.N2*(1-Bypass);
    block.Main.IC.O2 = Inlet.Air.O2*(1-Bypass) - 1.5*R(1) - .5*R(2) - .5*R(3);    
    
    flowName = fieldnames(Inlet.Air);
    for i = 1:1:length(flowName)
        if ~strcmp(flowName{i},'T')
            block.Bypass.IC.(flowName{i}) = Inlet.Air.(flowName{i})*Bypass;
        end
    end
    
    fuelName = fieldnames(Inlet.Fuel);
    for i = 1:1:length(fuelName) %% add other fuel species to main that are not in specinterest
        if ~strcmp(fuelName{i},'T') && nnz(strcmp(fuelName{i},specinterest))==0
            block.Main.IC.(fuelName{i}) = Inlet.Fuel.(fuelName{i});
        end
    end
    
    [~, H.air] = enthalpy(Inlet.Air);
    H.air1 = H.air*(1-Bypass);
    H.Bypass0 = H.air*Bypass;%sensible enthalpy at begining of bypass
    [~, H.fuel] = enthalpy(Inlet.Fuel);
    h = enthalpy(298);
    hrxn1 = 2*h.H2O+h.CO-h.CH4-1.5*h.O2; %Ch4 + 1.5 O2 --> CO + 2 H2O
    hrxn2 = h.CO2-h.CO-.5*h.O2; %CO + .5 O2 --> CO2 
    hrxn3 = h.H2O-h.H2-.5*h.O2; %H2 + .5 O2 -->  H2O
    H.Combust = hrxn1*R(1) + hrxn2*R(2) + hrxn3*R(3);
    
    [~, Y] = ode15s(@(t,y) SolveCombTemps(t,y,block,block.Main.IC,block.Bypass.IC,H), [0, 1e5], block.Scale(1:4));
    T = Y(end,:)';
    block.Pfactor = (NetFlow(block.Main.IC)+NetFlow(block.Bypass.IC))/block.Pdrop;
    %no dMdP or mFlow to update (fixed pressure drop)
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.Scale = [T;1;block.Pin.IC];
    block.IC = [ones(length(T),1);Bypass;1];%avoid divide by zero if bypass = 0
    block.Main.IC.T = T(1);
    block.Bypass.IC.T = T(2);
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Bypass = Y(5);
    Pin = Y(6);
    
    Inlet = checkSaturation(Inlet,block);
    NetOut = block.Pfactor*(Pin-Inlet.Pout);%total cold flow out
    spec = {'CH4','CO','CO2','H2','H2O'};
    for i =1:1:length(spec)
        if ~isfield(Inlet.Fuel,spec{i})
            Inlet.Fuel.(spec{i}) = 0;
        end
    end

    %merge flows
    inlets = {'Air';'Fuel'};
    H_in = 0;
    block.spec={};
    for j = 1:1:2
        H_in = H_in + enthalpy(Inlet.(inlets{j}));
        Spec = fieldnames(Inlet.(inlets{j}));
        Spec = Spec(~strcmp('T',Spec));
        if strcmp(inlets{j},'Air')
            b=Bypass;
        else b= 0;
        end
        for i = 1:1:length(Spec)
            if ismember(Spec{i},block.spec)
                ReactMix.(Spec{i}) = ReactMix.(Spec{i}) + Inlet.(inlets{j}).(Spec{i})*(1-b);
                if strcmp(inlets{j},'Air')
                    BypassFlow.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*b;
                end
            else
                block.spec(end+1) = Spec(i);
                ReactMix.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*(1-b);
                if strcmp(inlets{j},'Air')
                    BypassFlow.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*b;
                end
            end
        end
    end
    BypassFlow.T = Inlet.Air.T;

    %% 3 reaction:
    % CH4 + 1.5O2 --> CO + 2H2O
    % CO + .5O2 --> CO2
    % H2 + .5O2 --> H2O
    R.CH4 = ReactMix.CH4*block.Rxn1;
    R.CO = (ReactMix.CO+R.CH4)*block.Rxn2;
    R.H2 = ReactMix.H2*block.Rxn3;
    sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
    r = fieldnames(R);
    phi = sumR/(ReactMix.O2); %rich combustion limited by O2;
    EquivError = block.EquivSet - phi;
    if phi>1 %rich combustion
        for i = 1:1:length(r)
            R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
        end
    end

    for i = 1:1:length(block.spec)
        if strcmp(block.spec{i},'CH4')
            CombustMix.CH4 = ReactMix.CH4 - R.CH4;
        elseif strcmp(block.spec{i},'CO')
            CombustMix.CO = ReactMix.CO + R.CH4- R.CO;
        elseif strcmp(block.spec{i},'CO2')
            CombustMix.CO2 = ReactMix.CO2 + R.CO;
        elseif strcmp(block.spec{i},'H2')
            CombustMix.H2 = ReactMix.H2  - R.H2;
        elseif strcmp(block.spec{i},'H2O')
            CombustMix.H2O = ReactMix.H2O + 2*R.CH4 + R.H2;
        elseif strcmp(block.spec{i},'O2')
            CombustMix.O2 = ReactMix.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
        else
            CombustMix.(block.spec{i}) = ReactMix.(block.spec{i});
        end
    end
    CombustMix.T = Y(1);

    FlowFraction = NetOut/(NetFlow(CombustMix)+NetFlow(BypassFlow));%Outflow/Inflow
    Out.Main.T = Y(1);
    for i = 1:length(block.spec)
        Out.Main.(block.spec{i}) = CombustMix.(block.spec{i})*FlowFraction;%scale outflow by p-factor
    end

    spec = fieldnames(Inlet.Air);
    for i = 1:1:length(spec)
        Out.Bypass.(spec{i}) = BypassFlow.(spec{i})*FlowFraction;
    end
    Out.Bypass.T = Y(2);

    if strcmp(string1,'Outlet')
        Out.Pin = Pin;
        Tags.(block.name).EquivelanceRatio = phi;
        Tags.(block.name).Temperatures = Y(1:4)';
        Tags.(block.name).MassFlow = MassFlow(Out.Main)+MassFlow(Out.Bypass);
    elseif strcmp(string1,'dY')
        dY = 0*Y;

        Hcombusted = enthalpy(CombustMix);
        Hbypass = enthalpy(BypassFlow);

        Cp = SpecHeat(Out.Main);
        Cp2 = SpecHeat(Inlet.Air);

        Q_CombWall =(Y(1) - Y(3))*block.AmbConv*block.SurfA/1000;
        Q_BypassCasing = (Y(2) - Y(4))*block.SurfA*block.AmbConv/1000; %convection from bypass to casing
        Q_WallBypassC = (Y(3) - Y(2))*block.SurfA*block.AmbConv/1000;%Convection from Wall to Bypass
        Q_WallBypassR = ((Y(3))^4 - (Y(2))^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000;%Radiation from Wall to Bypass
        Q_CasingAmbR = ((Y(4))^4 - (block.Tamb)^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000; % Heat transfer from casing to ambient due to radiation
        Q_CasingAmbC = (Y(4) - block.Tamb)*block.AmbConv*block.SurfA/1000; % Heat transfer from casing to ambient due to convection

        dY(1) = (H_in - Hbypass - Hcombusted - Q_CombWall)*block.Ru*Y(1)/(Cp*block.Vol*Inlet.Pout); %main air flow through core combustor
        dY(2) = (Hbypass*FlowFraction - enthalpy(Out.Bypass) -Q_BypassCasing + Q_WallBypassR + Q_WallBypassC)*Y(2)*block.Ru/(Cp2*block.Vol*Inlet.Pout); %total change of bypass temp
        dY(3) =(Q_CombWall - Q_WallBypassR - Q_WallBypassC)/(block.MassWall * block.SpecHeatWall);%total change in wall temp
        dY(4) =(Q_BypassCasing - Q_CasingAmbR - Q_CasingAmbC)/(block.MassCasing * block.SpecHeatCasing); %Total change in Casing temp
        dY(5) = EquivError;
        if (Pin-Inlet.Pout)<0 && NetFlow(ReactMix)<=0
            dY(6) = -(Pin-Inlet.Pout);
        else
            dY(6) = (NetFlow(CombustMix)+NetFlow(BypassFlow)-NetOut)*block.Ru*Y(1)/(block.Vol);
        end
        Out = dY;
    end
end
end%Ends function Combustor

function dY = SolveCombTemps(t,Y,block,Main,Bypass,H)
dY = 0*Y;
Main.T = Y(1);
Bypass.T = Y(2);

Q_CombWall =(Y(1) - Y(3))*block.AmbConv*block.SurfA/1000;
Q_BypassCasing = (Y(2) - Y(4))*block.SurfA*block.AmbConv/1000; %convection from bypass to casing
Q_WallBypassC = (Y(3) - Y(2))*block.SurfA*block.AmbConv/1000;%Convection from Wall to Bypass
Q_WallBypassR = ((Y(3))^4 - (Y(2))^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000;%Radiation from Wall to Bypass
Q_CasingAmbR = ((Y(4))^4 - (block.Tamb)^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000; % Heat transfer from casing to ambient due to radiation
Q_CasingAmbC = (Y(4) - block.Tamb)*block.AmbConv*block.SurfA/1000; % Heat transfer from casing to ambient due to convection

[~,Hout] = enthalpy(Main);
[~, H_Bypass] = enthalpy(Bypass);%Sensible enthalpy at end of bypass

Cp = SpecHeat(Main);

dY(1) = (H.air1+H.fuel-H.Combust-Hout-Q_CombWall)/(Cp*NetFlow(Main));
dY(2) = (H.Bypass0 + Q_WallBypassR + Q_WallBypassC - H_Bypass - Q_BypassCasing)/(Cp*NetFlow(Bypass)); %total change of bypass temp
dY(3) = (Q_CombWall - Q_WallBypassR - Q_WallBypassC)/(block.MassWall*block.SpecHeatWall);%total change in wall temp
dY(4) = (Q_BypassCasing - Q_CasingAmbR - Q_CasingAmbC)/(block.MassCasing * block.SpecHeatCasing); %Total change in Casing temp  
end%Ends function SolveCombTemps