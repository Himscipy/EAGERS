function Out = Electrolyzer(varargin)
%This function models an electorlyzer it takes input Y states and returns dY: 
%Electrolyzer model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [reformer]), Cathode species ( [ CO2, H2O], N2, O2) anode species (CH4, CO, CO2, H2, H2O, N2, O2), [Reformer species] [Rate of internal reforming reactions] Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout'}
% Seven (7) outlets: {'Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureVoltage','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% for an elextrolyzer, Flow 1 is the cathode (steam), and Flow2 is the anode (cooling/heating air)
% Current is negative
global Tags
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.F=96485.339; % %Faraday's constant in Coulomb/mole
    block.Ru = 8.314472; % Universal gas constant in kJ/K*kmol
    block.nodes = block.rows*block.columns;
    %% Secondary Design Variables
    %%--Geometric Varibles  %
    block.t_plate1 = 0.003;                    % [m] thickness of plate1
    block.t_plate1_wall =0.005;                % [m] Thickness of the channel wall of the plate
    block.H_plate1channels = 0.002;            % [m] height of cathode channel
    block.W_plate1channels = 0.005;            % [m] width of channel                       
    block.H_plate2channels = 0.002;            % [m] height of anode channel
    block.W_plate2channels = 0.005;            % [m] width of  channel
    block.t_plate2 = 0.003;                    % [m] Thickness ofPlate
    block.t_plate2_wall=0.005;                 % [m] Thickness of the channel wall 
    block.Nu_flow1 = 4;                        %Nusselt # for square channel aspect ratio =3 & uniform temp
    block.Nu_flow2 = 4;                        %Nusselt # for square channel aspect ratio =3 & uniform temp

    %%Electrochemical parameters %%%
    % H2 +1/2O2 --> H2O (Nernst Eo)
    %SOFC Conductivity - Implemented Equation  sigma = A*e^(-deltaG/RT)/T
    switch block.FCtype
        case {'SOFC';'SOEC'}
            block.ElecConst = 2e3; %(K/ohm*m) Electrolyte Constant  %default SOFC  = 9e7
            block.deltaG = 8.0e3; %(kJ/kmol)
            block.t_Membrane = 18e-6;                     % [m] thickness of membrane
            block.t_Cath = 800e-6;                        % [m] thickness of cathode structure
            block.t_An = 50e-6;                           % [m] thickness of Anode structure
            block.t_Elec = block.t_Membrane+block.t_Cath+block.t_An;        % [m] thickness of complete electrolyte
    end

    %%Electrolyte
    block.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
    block.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
    block.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 

    %%--- Plate 1 -------%
    block.Density_plate1 = 2000;                                % [kg/m3]     density 
    block.C_plate1 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    block.k_plate1 = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    block.k_flow1 = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%-------Cathode Gas Stream---------%
    block.k_flow2 = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%---- Plate 2 -------%   
    block.Density_plate2 = 2000;                                % [kg/m3]     density 
    block.C_plate2 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    block.k_plate2 = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate


%%%---%%% end of user defined variables    
    %Dimensions
    block.A_Cell = block.L_Cell*block.W_Cell; %Cell Area
    block.A_Node = block.A_Cell/block.nodes; %node Area
    block.L_node = block.L_Cell/block.columns; %Node length in meters
    block.W_node = block.W_Cell/block.rows;  %Node Width in meters
    block.Dh_Flow1 = 4*(block.H_plate1channels*block.W_plate1channels)/(2*(block.H_plate1channels+block.W_plate1channels)); %(m) hydraulic diameter of channel
    block.Dh_Flow2 = 4*(block.H_plate2channels*block.W_plate2channels)/(2*(block.H_plate2channels+block.W_plate2channels)); %(m) hydraulic diameter of channel
    block.CH_Flow1 = block.W_node/(block.W_plate1channels+block.t_plate1_wall); %Number of channels in each node of the anode
    block.CH_Flow2 = block.W_node/(block.W_plate2channels+block.t_plate2_wall); %Number of channels in each node of the cathode 
    %% --- Plate 1 -------%
    block.A_plate1_elecCond = block.t_plate1_wall*block.L_node*block.CH_Flow1;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    block.A_plate1_heatCond = (block.H_plate1channels*block.t_plate1_wall + (block.W_plate1channels+block.t_plate1_wall)*block.t_plate1)*block.CH_Flow1; %[m^2] conduction area between nodes
    block.L_plate1_heatCond = block.H_plate1channels;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    block.Mass_plate1 = (block.H_plate1channels*block.t_plate1_wall + (block.W_plate1channels+block.t_plate1_wall)*block.t_plate1)*block.CH_Flow1*block.L_node*block.Density_plate1;
    %% ----- Cathode Gas Stream-----------%
    block.flow1_crossArea = block.H_plate1channels*block.W_plate1channels*block.CH_Flow1;              % [m2] Crossectional Area of Anode entrance
    block.h_flow1          = block.Nu_flow1*block.k_flow1/block.Dh_Flow1;                     % [W/m2/K]  Convection coefficient between the cathode gas and the Fuel Seperator plate
    block.A_flow1_plate1      = (2*block.H_plate1channels + block.W_plate1channels)*block.L_node*block.CH_Flow1;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    block.A_flow1_elec     = (block.W_plate1channels)*block.L_node*block.CH_Flow1;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    block.Vol_flow1        = block.H_plate1channels*block.W_plate1channels*block.L_node*block.CH_Flow1;               % [m3]  control volume
    %% --------Electrolyte-------------%
    block.A_Elec_Cond =  block.W_node*block.t_Elec;                   % [m2] Conduction surface area of electrolyte
    block.A_Elec_Heat_Cond = block.W_node*block.t_Elec;                    % [m2] Conduction surface area of electrolyte
    block.Vol_Elec = block.t_Elec*block.L_node*block.W_node;              % [m3] volume of electrolyte    
    %% ------- Anode Gas Stream---------%
    block.flow2_crossArea= block.H_plate2channels*block.W_plate2channels*block.CH_Flow2;       % [m2] Crossectional Area of Cathode entrance
    block.h_flow2= block.Nu_flow2*block.k_flow2/block.Dh_Flow2;                 % [W/m2/K]  Convection coefficient between the Anode gas and the Fuel Seperator plate
    block.A_flow2_plate2 = (2*block.H_plate2channels + block.W_plate2channels)*block.L_node*block.CH_Flow2;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    block.A_flow2_elec = block.W_plate2channels*block.CH_Flow2*block.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    block.Vol_flow2 = block.H_plate2channels*block.W_plate2channels*block.CH_Flow2*block.L_node;            % [m3]  control volume Cathode
    %% ---- Plate 2 -------%   
    block.A_plate2_elecCond = block.t_plate2_wall*block.L_node*block.CH_Flow2;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    block.A_plate2_heatCond = (block.H_plate2channels*block.t_plate2_wall + (block.W_plate2channels+block.t_plate2_wall)*block.t_plate2)*block.CH_Flow2; %[m^2] conduction area between nodes
    block.L_plate2_heatCond=block.H_plate2channels;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte	
    block.Mass_plate2 = (block.H_plate2channels*block.t_plate2_wall + (block.W_plate2channels+block.t_plate2_wall)*block.t_plate2)*block.L_node*block.CH_Flow2*block.Density_plate2;
    %% Pressure
    block.Flow1_Pout = block.PressureRatio*101;
    block.Flow1_Pinit = block.Flow1_Pout + block.Flow1Pdrop;
    block.Flow2_Pout =  block.PressureRatio*101;
    block.Flow2_Pinit = block.Flow2_Pout + block.Flow2Pdrop;
    
    switch block.Reformer %% Load mask parameters (flow direction)
        case 'methanator'
            block = FlowDir(block,3);
        case {'none'}
            block = FlowDir(block,2);
    end
    %estimate # of cells
    if block.ClosedCathode
        block.Cells = ceil(block.RatedStack_kW*1000/(1.3*1e4*block.L_Cell*block.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
    else
        if strcmp(block.Specification,'cells')
            block.Cells = block.SpecificationValue;
            block.Specification = 'power density';
            block.SpecificationValue = block.RatedStack_kW*100/(block.L_Cell*block.W_Cell*block.Cells);
        elseif strcmp(block.Specification,'power density')
            block.Cells = ceil(block.RatedStack_kW*100/(block.L_Cell*block.W_Cell*block.SpecificationValue)); %# of cells in stack
        elseif strcmp(block.Specification,'current density')
            block.Cells = ceil(block.RatedStack_kW*1000/(1.3*1e4*block.L_Cell*block.W_Cell*block.SpecificationValue)); %# of cells in stack (assumes voltage of 1.3)
        elseif strcmp(block.Specification,'voltage')
            block.Cells = ceil(block.RatedStack_kW*1000/(block.SpecificationValue*1e4*block.L_Cell*block.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
        end 
    end
     %% Estimate heat generated
    [h,~] = enthalpy(block.TpenAvg,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    Vbalance = 1./(2*block.F).*h_rxn3; %voltage that balances heat
    if block.ClosedCathode
        block.SpecificationValue = Vbalance;
        block.Voltage = block.SpecificationValue;
        block.Specification = 'voltage';
    else
        if strcmp(block.Specification,'power density')
            block.Voltage = 1.3;
        elseif strcmp(block.Specification,'current density')
            block.Voltage = block.RatedStack_kW/block.Cells*1000/(block.A_Cell*(100^2))/block.SpecificationValue; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        elseif strcmp(block.Specification,'voltage')
            block.Voltage = block.SpecificationValue;
        end
    end
    
    %% %% 1st guess at Initial Condition
    Current = zeros(block.nodes,1);
    block.Current = zeros(block.nodes,1);
    if strcmp(block.Specification,'power density')
        i_avg = -block.SpecificationValue/block.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(block.Specification,'current density')
        i_avg = -block.SpecificationValue;
    elseif strcmp(block.Specification,'voltage')
        i_avg = -block.RatedStack_kW/block.Cells*1000/(block.A_Cell*(100^2))/block.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:block.rows
        Current(1+block.columns*(j-1):block.columns*j) =linspace(2,1,block.columns)/sum(linspace(2,1,block.columns))*i_avg*(100^2)*block.A_Cell/block.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end
    block.Current.CO = 0*Current;
    block.Current.H2 = Current;
    
    block.StackTempIn = block.TpenAvg-.9*block.deltaTStack;
    block.T.Flow1 = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Elec = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Flow2 = zeros(block.nodes,1) + block.TpenAvg;
   

    Flow2In = block.Flow2Spec;
    Flow2In.T = block.TpenAvg;
    Cp = SpecHeat(Flow2In);
    Q = (block.Voltage-Vbalance)*sum(abs(block.Current.H2 + block.Current.CO))*block.Cells/1000; %estimated kW of heat generated in stack
    block.AirFlow = abs(Q)/(Cp*block.deltaTStack); %estimate of cooling air flow

    block.F2Spec = unique([{'O2';};fieldnames(block.Flow2Spec)]);
    
    block.SteamFlow  = sum(abs(Current))/(2*block.F*1000)/(block.H2O_Utilization*block.Flow1Spec.H2O)*block.Cells; % H2O flow rate,  current/(2*block.F*1000) = kmol H2
    SteamSpec = fieldnames(block.Flow1Spec);
    block.F1Spec = unique([{'H2';'H2O';};SteamSpec]);
    
    Inlet = InletFlow(block);
    if strcmp(block.Specification,'current density')
        Inlet.NetCurrent = i_avg*(block.A_Cell*(100^2));
    end
    Inlet.Flow1.T = block.TpenAvg -.75*block.deltaTStack;
    Inlet.Flow2.T = block.TpenAvg -.75*sign(Q)*block.deltaTStack;
    
    %% Run Initial Condition
    [Flow1, Flow2,block,Inlet] = solveInitCond(Inlet,block,1);
    
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout'};
    block.NetCurrent.IC = sum(block.Current.H2 + block.Current.CO);
    block.Flow1.IC = Inlet.Flow1;
    block.Flow2.IC =  Inlet.Flow2;
    block.Flow1Pout.IC = block.Flow1_Pout;
    block.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.Flow2Pout.IC = block.Flow2_Pout;
    block.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTflow1','MeasureTflow2'};
    block.Flow1Out.IC  = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    block.Flow2Out.IC =  MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.MeasureVoltage.IC = block.Voltage;
    block.MeasurePower.IC = sum(abs(block.Current.H2 + block.Current.CO)*block.Voltage*block.Cells)/1000;%power in kW
    block.MeasureTpen.IC = block.T.Elec;
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));

    block.P_Difference = {'Flow1Pin','Flow1Pout'; 'Flow2Pin', 'Flow2Pout';};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    if strcmp(block.Specification,'power density')
        block.Specification = 'current density'; %convergence of power done by controller initialization
    end
    
    Flow1SpecNew = fieldnames(Inlet.Flow1);
    Flow1SpecAll = unique([block.F1Spec;Flow1SpecNew]);
    Flow1SpecAll = Flow1SpecAll(~strcmp('T',Flow1SpecAll));
    for i = 1:1:length(Flow1SpecAll)
        if ~ismember(Flow1SpecAll{i},Flow1SpecNew)
            Inlet.Flow1.(Flow1SpecAll{i})=0;
        end
    end
    block.F1Spec = Flow1SpecAll;
    
    Flow2SpecNew = fieldnames(Inlet.Flow2);
    Flow2SpecAll = unique([block.F2Spec;Flow2SpecNew]);
    Flow2SpecAll = Flow2SpecAll(~strcmp('T',Flow2SpecAll));
    for i = 1:1:length(Flow2SpecAll)
        if ~ismember(Flow2SpecAll{i},Flow2SpecNew)
            Inlet.Flow2.(Flow2SpecAll{i})=0;
        end
    end
    block.F2Spec = Flow2SpecAll;

    block.Flow1_Pinit = Inlet.Flow1Pout + block.Flow1Pdrop;
    block.Flow2_Pinit = Inlet.Flow2Pout + block.Flow2Pdrop;
    block.Flow1Pout.IC = Inlet.Flow1Pout;
    block.Flow2Pout.IC = Inlet.Flow2Pout;
    block.StackPower = abs(Inlet.NetCurrent)*block.Cells/1000*block.Voltage;
    %%--%%
    [Flow1, Flow2,block,~] = solveInitCond(Inlet,block,2);
    %%%
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one
    block.Flow1Out.IC  = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    block.Flow2Out.IC = MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.MeasureVoltage.IC = block.Voltage;
    block.MeasurePower.IC = sum(abs(block.Current.H2 + block.Current.CO)*block.Voltage*block.Cells)/1000;%power in kW
    block.MeasureTpen.IC = block.T.Elec;
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));
    Out = block;
else%running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    nodes = block.nodes;
    %% add species that may not be in inlet
    inFields = fieldnames(Inlet.Flow1);
    for i = 1:1:length(block.F1Spec)
        if ~ismember(block.F1Spec{i},inFields)
            Inlet.Flow1.(block.F1Spec{i}) = 0;
        end
    end
    inFields = fieldnames(Inlet.Flow2);
    for i = 1:1:length(block.F2Spec)
        if ~ismember(block.F2Spec{i},inFields)
            Inlet.Flow2.(block.F2Spec{i}) = 0;
        end
    end
    %% separate out temperatures
    T_plate1 = Y(1:nodes);
    Flow1.Outlet.T = Y(nodes+1:2*nodes);
    T_Elec = Y(2*nodes+1:3*nodes);
    Flow2.Outlet.T = Y(3*nodes+1:4*nodes); 
    T_plate2 = Y(4*nodes+1:5*nodes); 
    switch block.Reformer
        case 'methanator'
            nT = 6*nodes; % # of temperature states
            Flow3.Outlet.T = Y(5*nodes+1:6*nodes); %only gets used if internal Methanator exists, otherwise these values are actually QT1
        case 'none'
            nT = 5*nodes; % # of temperature states
    end
    n = nT;

    %Current
    nCurrent = Y(end-nodes-1:end-2);
    %Pressure
    P_flow1 = Y(end-1); %pressure
    P_flow2 = Y(end); %pressure

    %% Cathode
    for i = 1:1:length(block.F1Spec)
        Flow1.Outlet.(block.F1Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
    end
    for j = 1:1:length(block.Flow1Dir(1,:));
        if j==1%first column recieves fresh inlet
            k = block.Flow1Dir(:,1);
            Flow1.Inlet.T(k,1) = Inlet.Flow1.T; 
            for i = 1:1:length(block.F1Spec)
                Flow1.Inlet.(block.F1Spec{i})(k,1) = Inlet.Flow1.(block.F1Spec{i})/block.Cells/length(k); 
            end
        else%subsequent columns recieve outlet of previous column
            k2 = block.Flow1Dir(:,j);
            Flow1.Inlet.T(k2,1) = Flow1.Outlet.T(k);
            for i = 1:1:length(block.F1Spec)
                Flow1.Inlet.(block.F1Spec{i})(k2,1) = Flow1.Outlet.(block.F1Spec{i})(k); 
            end
            k = k2;
        end
    end
    Flow1Out.T  = mean(Flow1.Outlet.T(block.Flow1Dir(:,end))); %temperature 
    for i = 1:1:length(block.F1Spec)
        Flow1Out.(block.F1Spec{i}) = max(0,sum(Flow1.Outlet.(block.F1Spec{i})(block.Flow1Dir(:,end)))*block.Cells);%avoid sending negative outlets
    end

    %% Anode 
    for i = 1:1:length(block.F2Spec)
        Flow2.Outlet.(block.F2Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
    end
    for j = 1:1:length(block.Flow2Dir(1,:))
        k2 =block.Flow2Dir(:,j);
        if j==1 % first column of fuel flow direction
            Flow2.Inlet.T(k2,1) = Inlet.Flow2.T;
            for i = 1:1:length(block.F2Spec)
                Flow2.Inlet.(block.F2Spec{i})(k2,1) = Inlet.Flow2.(block.F2Spec{i})/block.Cells/length(k2);
            end
        else
            Flow2.Inlet.T(k2,1) = Flow2.Outlet.T(k);
            for i = 1:1:length(block.F2Spec)
                Flow2.Inlet.(block.F2Spec{i})(k2,1) = Flow2.Outlet.(block.F2Spec{i})(k);
            end
        end
        k = k2;
    end
    Flow2Out.T  = mean(Flow2.Outlet.T(block.Flow2Dir(:,end))); %temperature 
    for i = 1:1:length(block.F2Spec)
        Flow2Out.(block.F2Spec{i}) = max(0,sum(Flow2.Outlet.(block.F2Spec{i})(block.Flow2Dir(:,end)))*block.Cells);%avoid sending negative outlets
    end

    %% Methanator
    switch block.Reformer
        case 'methanator' %secondary inlet introduces CO2 stream
            for i = 1:1:length(block.F2Spec)
                Flow3.Outlet.(block.F2Spec{i}) = Y(n+1:n+nodes); n = n+nodes;
            end
            for j = 1:1:length(block.Flow3Dir(1,:))
                if j==1
                    k = block.Flow3Dir(:,1);
                    Flow3.Inlet.T(k,1) = Inlet.CarbonDioxide.T;
                    for i = 1:1:length(block.F2Spec)
                        Flow3.Inlet.(block.F2Spec{i})(k,1) = Flow2.Outlet.(block.F2Spec{i})(block.Flow2Dir(:,1))*block.MethSpacing;
                        if isfield(Inlet.CarbonDioxide,block.F2Spec{i})
                            Flow3.Inlet.(block.F2Spec{i})(k,1) = Flow3.Inlet.(block.F2Spec{i})(k,1) +Inlet.CarbonDioxide.(block.F2Spec{i})/block.Cells/length(k)*block.MethSpacing;
                        end
                    end
                else
                    k2 = block.Flow3Dir(:,j);
                    Flow3.Inlet.T(k2,1) = Flow3.Outlet.T(k);
                    for i = 1:1:length(block.F2Spec)
                        Flow3.Inlet.(block.F2Spec{i})(k2,1) = Flow3.Outlet.(block.F2Spec{i})(k);
                    end
                    k = k2;
                end
            end
    end
    %%Nernst & Losses
    FuelCellNernst(Flow1,Flow2,nCurrent,T_Elec,P_flow2,block)
    Voltage =  sum(Tags.(block.name).nVoltage'.*(nCurrent/sum(nCurrent)));
    Current.CO = Tags.(block.name).I_CO';
    Current.H2 = Tags.(block.name).I_H2';
    [h,~] = enthalpy(T_Elec,{'H2','H2O','O2','CO','CO2'});
    nPower = Voltage*nCurrent/1000; %cell power in kW
    Qgen = Current.H2/(2000*block.F).*(h.H2+.5*h.O2-h.H2O) + Current.CO/(2000*block.F).*(h.CO+.5*h.O2-h.CO2) - nPower;%kW of heat generated by electrochemistry (per node & per cell)
    Power = Voltage*abs(Inlet.NetCurrent)*block.Cells/1000;
    if strcmp(string1,'Outlet')
        H2O_in = Inlet.Flow1.H2O;
        H2O_out = Flow1Out.H2O;
        %%Outlet Ports
        Out.Flow1Out  = Flow1Out;
        Out.Flow2Out = Flow2Out;
        Out.Flow1Pin = P_flow1;
        Out.Flow2Pin = P_flow2;
        Out.MeasureVoltage = Voltage;
        Out.MeasurePower = Power;
        Out.MeasureTpen = Y(2*nodes+1:3*nodes);
        Out.MeasureTflow1 = Y(nodes+block.Flow1Dir(:,end));
        Out.MeasureTflow2 = Y(3*nodes+block.Flow2Dir(:,end));
        %% Tags
        Tags.(block.name).Voltage = Voltage;
        Tags.(block.name).H2Outilization = (H2O_in - H2O_out)./H2O_in;
        Tags.(block.name).Tpen = T_Elec';
        Tags.(block.name).TcathOut = Flow1Out.T;
        Tags.(block.name).TanodeOut = Flow2Out.T;
        Tags.(block.name).Current = sum(abs(nCurrent));
        Tags.(block.name).StackdeltaT = Flow2Out.T-Inlet.Flow2.T;
        Tags.(block.name).PENavgT = sum(T_Elec)/block.nodes;
        Tags.(block.name).MaxPEN = max(T_Elec);
        Tags.(block.name).PENdeltaT = Tags.(block.name).MaxPEN-min(T_Elec);
        Tags.(block.name).dTdX = (T_Elec-T_Elec(block.HTadjacent(:,2)))'/(block.L_Cell/block.columns);
        Tags.(block.name).dTdY = (T_Elec-T_Elec(block.HTadjacent(:,4)))'/(block.W_Cell/block.rows);
        Tags.(block.name).MaxdTdX = max(abs([Tags.(block.name).dTdX;Tags.(block.name).dTdY;]));
        Tags.(block.name).Q_gen = sum(Qgen*block.Cells); %kW of heat generated by electrochemistry
        Tags.(block.name).Efficiency = (NetFlow(Flow1Out)*HeatingValue(Flow1Out))/Power;
    elseif strcmp(string1,'dY')  
%         Voltage =  Tags.(block.name).Voltage;
        switch block.Reformer
            case 'methanator'
                MethCurrent.H2 = zeros(nodes,1);
                MethCurrent.CO = zeros(nodes,1);
                [Rref,RefOut] = KineticReformation(block.method,Flow3,P_flow1,block.KineticCoeff3,MethCurrent,block);%% Kinetic reaction rates  (WGS is always near equilibrium)
                [~,RefOut] = KineticReformation(block.method,Flow3,P_flow2,block.KineticCoeff3,zeros(block.nodes,1),block);%% Kinetic reaction rates  (WGS is always near equilibrium)
    %             R.CH4ref = Rref.CH4;
    %             R.WGSref = Rref.WGS;
        end

        
        switch block.FCtype%ion transport across membrane (total enthalpy)
            case {'SOFC';'SOEC'}
                Qion = nCurrent/(4000*block.F).*h.O2; %O2 ion crossing over (kW)
            case {'MCFC';'MCEC'}
                Qion = nCurrent/(4000*block.F).*h.O2 + nCurrent/(2000*block.F).*h.CO2;% O2 & CO2 ion crossing over
        end

        %% Q %% Heat transfer & Generation
        switch block.Reformer
            case 'methanator'
                QT = block.HTcond*Y(1:6*nodes) + block.HTconv*Y(1:6*nodes);
            case {'none'}
                QT = block.HTcond*Y(1:5*nodes) + block.HTconv*Y(1:5*nodes);
        end
        
        %energy flows & sepcific heats
        Hout2 = enthalpy(Flow2.Outlet);
        Hin2 = enthalpy(Flow2.Inlet);
        Hout1 = enthalpy(Flow1.Outlet);
        Hin1 = enthalpy(Flow1.Inlet);

        %% %% solve for dY in order of states
        dY = 0*Y;
        %%Temperatures
        dY(1:nodes)= QT(1:nodes)./block.tC(1:nodes);  %Cathode Plate
        for i=1:1:length(block.Flow1Dir(1,:)) %having the downstream nodes change temperature with the upstream nodes prevents propogation issues when taking larger time steps
            k = block.Flow1Dir(:,i);
            dY(nodes+k)= (QT(nodes+k) + Hin2(k) - Hout2(k) - Qion(k))./block.tC(nodes+k); 
            if i>1
                dY(nodes+k) = dY(nodes+k)+dY(nodes+kprev);
            end
            kprev = k;
        end
        dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./block.tC(2*nodes+1:3*nodes); %Electrolyte Plate
        for i=1:1:length(block.Flow2Dir(1,:))
            k = block.Flow2Dir(:,i);
            dY(3*nodes+k)= (QT(3*nodes+k) + Hin1(k) - Hout1(k) + Qion(k) - nPower(k) - Qgen(k))./block.tC(3*nodes+k);
            if i>1
                dY(3*nodes+k) = dY(3*nodes+k)+dY(3*nodes+kprev);
            end
            kprev = k;
        end
        dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./block.tC(4*nodes+1:5*nodes);

        n =nT;
        %%Cathode Species
        for i = 1:1:length(block.F1Spec)
            if strcmp(block.F1Spec{i},'H2O')
                dY(n+1:n+nodes)= (Flow1.Inlet.H2O - Flow1.Outlet.H2O + nCurrent/(2*block.F*1000))./block.tC(n+1:n+nodes);  %H2O species concentration with CO2 crossover
            elseif strcmp(block.F1Spec{i},'H2') 
                dY(n+1:n+nodes)= (Flow1.Inlet.H2 - Flow1.Outlet.H2 - nCurrent/(2*block.F*1000))./block.tC(n+1:n+nodes);%H2 species concentration with O2 crossover
            else
                dY(n+1:n+nodes)= (Flow1.Inlet.(block.F1Spec{i}) - Flow1.Outlet.(block.F1Spec{i}))./block.tC(n+1:n+nodes);%all other species concentration
            end
            n = n+nodes;
        end

        %%Anode Species
        for i = 1:1:length(block.F2Spec)
            if strcmp(block.F2Spec{i},'O2')
                dY(n+1:n+nodes)= (Flow2.Inlet.O2 - Flow2.Outlet.O2 - nCurrent/(4*block.F*1000))./block.tC(n+1:n+nodes);%O2 species concentration with O2 crossover
            else
                dY(n+1:n+nodes)= (Flow2.Inlet.(block.F2Spec{i}) - Flow2.Outlet.(block.F2Spec{i}))./block.tC(n+1:n+nodes); %all species concentration
            end
            n = n+nodes; 
        end

        %%Methanator
        switch block.Reformer
            case 'methanator'
                Hout3 = enthalpy(Flow3.Outlet);
                Hin3 = enthalpy(Flow3.Inlet);
                for i=1:1:length(block.ReformFlowDir(1,:))
                    k = block.ReformFlowDir(:,i);
                    dY(5*nodes+k)= (block.RefSpacing*QT(5*nodes+k) + Hin3(k) - Hout3(k))./block.tC(5*nodes+k);  %Fuel Methanator Channels
                    if i>1
                        dY(5*nodes+k) = dY(5*nodes+k)+dY(5*nodes+kprev);
                    end
                    kprev = k;
                end
                for i = 1:1:length(block.F2Spec)
                    dY(n+1:n+nodes)= (RefOut.(block.F2Spec{i}) - Flow3.Outlet.(block.F2Spec{i}))./block.tC(n+1:n+nodes);   %all species concentrations
                    n = n+nodes;
                end
        end
        %%Current % note sign convention of current is negative for electrolyzers !!
        dY(n+1:n+nodes) = (Inlet.NetCurrent - sum(nCurrent) + (Tags.(block.name).nVoltage' - Voltage)./Tags.(block.name).ASR.*(block.A_Cell*100^2))./block.tC(end-nodes-1:end-2); n = n+nodes; %error in A/cm^2 * area  
        
        

        %%Pressure
        Nanode = block.Pfactor2*max(0.01,(P_flow2-Inlet.Flow2Pout));%total anode flow out
        Ncath = block.Pfactor1*max(0.1,(P_flow1-Inlet.Flow1Pout));%total cathode flow out
        dY(n+1) = (NetFlow(Inlet.Flow1)-Ncath)*block.Ru*Inlet.Flow1.T/block.tC(n+1);%working with total flow rates so must multiply by nodes & cells
        dY(n+2) = (NetFlow(Inlet.Flow2)-Nanode)*block.Ru*Inlet.Flow2.T/block.tC(n+2);
        Out = dY;
    end
end
end%Ends function Electrolyzer

function [Flow1, Flow2,block,Inlet] = solveInitCond(Inlet,block,firstSolve)
%% loop to converge voltage 
% SOEC & MCEC, adjust current and H2 flow to achieve power density and H2O utilization
Tol = 1e-3;
error = 10*Tol;
count = 1;
Flow3 = []; %Currently not set up for methane production internally
block.R_CH4 = zeros(block.nodes,1);
block.R_WGS = zeros(block.nodes,1);
while abs(error)>Tol %iterate to find the initial current (holding power or voltage fixed)
    Flow1 = FCin2Out(block.T.Flow1,Inlet.Flow1,block.Flow1Dir, block.FCtype,block.Cells,block.Current,[],'cathode');
    Flow2 = FCin2Out(block.T.Flow2,Inlet.Flow2,block.Flow2Dir, block.FCtype,block.Cells,block.Current,[],'anode');
%     Utilization = (sum(abs(block.Current.H2 + block.Current.CO))*block.Cells/(2000*F))/(Inlet.Flow1.H2O);
    if count==1 && firstSolve==1
        [block.Tstates,block.HTcond,block.HTconv]= SteadyTemps(block,Inlet.Flow1,Inlet.Flow2);
    else
        [~, Y] = ode15s(@(t,y) DynamicTemps(t,y,block,Flow1,Flow2,Flow3), [0, 1e5], block.Tstates);
        block.Tstates = Y(end,:)';
    end
    %organize temperatures
    
    block.T.Elec =  block.Tstates(2*block.nodes+1:3*block.nodes);
    Tcorrection = block.TpenAvg - mean(block.T.Elec);
    block.T.Flow1 =  block.Tstates(1*block.nodes+1:2*block.nodes)+Tcorrection;
    block.T.Flow2 = block.Tstates(3*block.nodes+1:4*block.nodes)+Tcorrection;

    %% Nernst & Losses
    normTemp = block.T.Elec+Tcorrection; %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence)
    FuelCellNernst(Flow1,Flow2,block.Current,normTemp,block.Flow2_Pinit,block)
    
    %%Converge current distribution
    [block,error,scale] = redistributeCurrent(block,Inlet,count,firstSolve);
    TotCurrent = abs(sum(block.Current.H2 + block.Current.CO));
    if firstSolve ==1
        if strcmp(block.Specification,'voltage') || strcmp(block.Specification,'current density')
            block.Cells = ceil((block.RatedStack_kW*1000/block.Voltage)/(sum(-block.Current))); %re-calculate the # of cells
        end
        block.SteamFlow  = TotCurrent/(2*block.F*1000)/(block.H2O_Utilization*block.Flow1Spec.H2O)*block.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2

        if block.ClosedCathode
            block.AirFlow = 0;
            Inlet = InletFlow(block);
            Inlet.Flow2.T = block.TpenAvg; % no anode flow in
        else
            Q_cathode = block.Cells*sum(NetFlow(Flow1.Outlet).*SpecHeat(Flow1.Outlet).*(Flow1.Outlet.T - Flow1.Inlet.T));
            [h,~] = enthalpy(block.TpenAvg,{'H2','H2O','O2'});
            h_rxn3 = h.H2+.5*h.O2-h.H2O;
            Vbalance = 1/(2*block.F)*h_rxn3; %voltage that balances heat
            block.AirFlow = abs((block.Cells*(block.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)/(33*50); %air flow is extra heat / Cp* deltaT
            Inlet = InletFlow(block);
            if ((block.Cells*(block.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)>0
                Inlet.Flow2.T = block.TpenAvg-100; %cooling stack
            else
                Inlet.Flow2.T= block.TpenAvg+100;%heating stack
            end
        end
        Inlet.Flow1.T = block.TpenAvg -.75*block.deltaTStack;
        
    end
    count= count+1;
end
%% Finish initialization by organizing state variables
block.Pfactor2 = block.AirFlow/block.Flow2Pdrop;
block.Pfactor1 = block.SteamFlow/block.Flow1Pdrop;
block = Set_IC(block,Flow1,Flow2,Flow3);
end%Ends function SolveInitCond

function Inlet = InletFlow(block) %only used 1st time through initialization (before we know what is connected to inlet
% Anode (oxygen production)
if block.AirFlow>0
    for i = 1:1:length(block.F2Spec)
        if isfield(block.Flow2Spec,block.F2Spec{i})
            Inlet.Flow2.(block.F2Spec{i}) = block.Flow2Spec.(block.F2Spec{i})*block.AirFlow;%flow rate of every species entering the anode (or reformer if there is one)
        else Inlet.Flow2.(block.F2Spec{i}) = 0;
        end
    end
else % no dilution air
    for i = 1:1:length(block.F2Spec)
        Inlet.Flow2.(block.F2Spec{i}) = 0;
    end
end

%Cathode H2 production
switch block.FCtype
    case {'MCEC';'SOEC'}
        for i = 1:1:length(block.F1Spec)
            if isfield(block.Flow1Spec,block.F1Spec{i})
                Inlet.Flow1.(block.F1Spec{i}) = block.Flow1Spec.(block.F1Spec{i})*block.SteamFlow;
            else Inlet.Flow1.(block.F1Spec{i}) = 0;
            end
        end
end
end%Ends function InletFlow

function block = Set_IC(block,Flow1,Flow2,Flow3)
Cp_2 = SpecHeat(Flow2.Inlet);
Cp_1 = SpecHeat(Flow1.Outlet);
if ~isempty(Flow3)
    Cp.meth = SpecHeat(Flow3.Outlet);
end
switch block.Reformer
    case 'none'
        NumOfStates = (5 + length(block.F2Spec) + length(block.F1Spec) + 1)*block.nodes + 2; % 5 temperatures, anode species & cathode species & current at each node and 2 states for anode/cathode pressure 
    case 'methanator'
end

block.Scale = ones(NumOfStates,1); %2 states for anode/cathode pressure
block.IC = block.Scale;
block.UpperBound = inf*ones(NumOfStates,1);
block.LowerBound = zeros(NumOfStates,1); %need to make this -inf for current states 
block.tC = block.IC; % time constant for derivative dY
block.Scale(1:5*block.nodes) = block.Tstates(1:5*block.nodes);%temperature (K)

block.tC(1:block.nodes) = (block.Mass_plate2*block.C_plate2);
block.tC(1+block.nodes:2*block.nodes) = (block.Vol_flow1*Cp_2*block.Flow1_Pinit./(block.Ru*block.T.Flow1));
block.tC(2*block.nodes+1:3*block.nodes) = (block.Vol_Elec*block.Density_Elec*block.C_Elec);
block.tC(3*block.nodes+1:4*block.nodes) = (block.Vol_flow2*Cp_1*block.Flow2_Pinit./(block.Ru*block.T.Flow2));
block.tC(4*block.nodes+1:5*block.nodes) = (block.Mass_plate1*block.C_plate1);

n = 5*block.nodes;
for i = 1:1:length(block.F1Spec)
    block.tC(n+1:n+block.nodes) = (block.Vol_flow1*block.Flow1_Pinit)./(block.T.Flow1*block.Ru);  % cathode 
    if any(Flow1.Outlet.(block.F1Spec{i})==0)
        block.IC(n+1:n+block.nodes) = Flow1.Outlet.(block.F1Spec{i})./NetFlow(Flow1.Outlet);
        block.Scale(n+1:n+block.nodes) = NetFlow(Flow1.Outlet); n = n+block.nodes; %cathode flows
    else
        block.Scale(n+1:n+block.nodes) = Flow1.Outlet.(block.F1Spec{i}); n = n+block.nodes; %cathode flows
    end
end

for i = 1:1:length(block.F2Spec)
    block.tC(n+1:n+block.nodes) = (block.Vol_flow2*block.Flow2_Pinit)./(block.T.Flow2*block.Ru); %anode
    if any(Flow2.Outlet.(block.F2Spec{i})==0) %concentration less than 1%
        block.IC(n+1:n+block.nodes) = Flow2.Outlet.(block.F2Spec{i})./NetFlow(Flow2.Outlet);%concentration
        block.Scale(n+1:n+block.nodes) = NetFlow(Flow2.Outlet); %anode flow
    else
        block.Scale(n+1:n+block.nodes) = Flow2.Outlet.(block.F2Spec{i}); %individual species flow
    end   
    n = n+block.nodes;
end

switch block.Reformer
    case 'methanator'
        for i = 1:1:length(block.F1Spec)
            block.tC(n+1:n+block.nodes) = (block.Vol_flow3*block.FuelPinit)./(block.T.Flow3*block.Ru);
            if any(Flow3.Outlet.(block.F1Spec{i})==0)
                block.IC(n+1:n+block.nodes) = Flow3.Outlet.(block.F1Spec{i})./NetFlow(Flow3.Outlet);
                block.Scale(n+1:n+block.nodes) = NetFlow(Flow3.Outlet); 
            else
                block.Scale(n+1:n+block.nodes) = Flow3.Outlet.(block.F1Spec{i}); 
            end
            n = n+block.nodes;
        end
end
% note sign convention of current is negative for electrolyzers !!
block.IC(n+1:n+block.nodes) = -1;  %current
block.tC(n+1:n+block.nodes) = 1;  %current
block.Scale(n+1:n+block.nodes) = abs(block.Current.H2+block.Current.CO); %current 
block.LowerBound(n+1:n+block.nodes) = -inf; n = n+block.nodes; %current

block.tC(n+1) = (block.Vol_flow2*block.nodes*block.Cells);  %pressure
block.tC(n+2) = (block.Vol_flow1*block.nodes*block.Cells); %pressure
block.Scale(n+1) = block.Flow1_Pinit;%pressure
block.Scale(n+2) = block.Flow2_Pinit;%pressure
end%Ends function Set_IC

function dY = DynamicTemps(t,Y,block,Flow1,Flow2,Flow3)
dY = 0*Y;
nodes = block.nodes;

if isfield(block,'tC')
    tC = block.tC(1:5*nodes);
else
    Cp_1 = 42;% kJ/kmol*K
    Cp_2 = 33;% kJ/kmol*K
    tC(1:nodes,1) = (block.Mass_plate1*block.C_plate1);
    tC(1+nodes:2*nodes,1) = (block.Vol_flow1*Cp_1*block.Flow1_Pinit./(block.Ru*block.T.Flow1));
    tC(2*nodes+1:3*nodes,1) = (block.Vol_Elec*block.Density_Elec*block.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (block.Vol_flow2*Cp_2*block.Flow2_Pinit./(block.Ru*block.T.Flow2));
    tC(4*nodes+1:5*nodes,1) = (block.Mass_plate2*block.C_plate2);
end

h = enthalpy(Y(1+2*nodes:3*nodes),{'H2','H2O','O2','CO','CO2'});
Power = block.Voltage*(block.Current.H2+block.Current.CO)/1000; %cell power in kW
Qgen = block.Current.H2/(2000*block.F).*(h.H2+.5*h.O2-h.H2O) + block.Current.CO/(2000*block.F).*(h.CO+.5*h.O2-h.CO2)-Power;%kW of heat generated by electrochemistry (per node & per cell)
switch block.FCtype%ion transport across membrane (total enthalpy)
    case {'SOFC';'SOEC'}
        Qion = (block.Current.H2+block.Current.CO)/(4000*block.F).*h.O2; %O2 ion crossing over (kW)
    case {'MCFC';'MCEC'}
        Qion = (block.Current.H2+block.Current.CO)*(1/(4000*block.F).*h.O2 + 1/(2000*block.F).*h.CO2);% O2 & CO2 ion crossing over
end

QT = block.HTcond*Y + block.HTconv*Y;

Flow2.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    if j~=1
        Flow2.Inlet.T(k,1) = Flow2.Outlet.T(kprev);
    end
    kprev = k;
end

Flow1.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
    k = block.Flow2Dir(:,j);
    if j~=1
        Flow1.Inlet.T(k,1) = Flow1.Outlet.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
Hout1 = enthalpy(Flow1.Outlet);
Hin1 = enthalpy(Flow1.Inlet);
Hout2 = enthalpy(Flow2.Outlet);
Hin2 = enthalpy(Flow2.Inlet);

switch block.Reformer
    case 'methanator'
        Flow3.Outlet.T = Y(5*nodes+1:6*nodes);
        for j = 1:1:length(block.Flow3Dir(1,:));%1:columns
            k = block.Flow3Dir(:,j);
            if j~=1
                Flow3.Inlet.T(k,1) = Flow3.Outlet.T(kprev);
            end
            kprev = k;
        end
        Hout3 = enthalpy(Flow3.Outlet);
        Hin3 = enthalpy(Flow3.Inlet);
end

if block.ClosedCathode %%energy balance
    Qimbalance = sum((Hin2(block.Flow1Dir(:,1))) - sum(Hout2(block.Flow1Dir(:,end)))) + sum(Hin1(block.Flow2Dir(:,1)))  - sum(Hout1(block.Flow2Dir(:,end))) - sum(Power);
    Power = Power + Qimbalance*Power./sum(Power);
    Qgen = block.Current/(2000*block.F).*(h.H2+.5*h.O2-h.H2O) - Power;%kW of heat generated by electrochemistry (per node & per cell)
end

dY(1:nodes)= QT(1:nodes)./tC(1:nodes);  %Ox Sep Plate
dY(1+nodes:2*nodes)= (QT(1+nodes:2*nodes) + Hin2 - Hout2 - Qion)./tC(1+nodes:2*nodes); %Cathode
dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./tC(2*nodes+1:3*nodes); %Electrolyte Plate
dY(1+3*nodes:4*nodes)= (QT(1+3*nodes:4*nodes) + Hin1 - Hout1 + Qion - Power - Qgen)./tC(1+3*nodes:4*nodes);  %Anode
dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
switch block.Reformer
    case 'methanator'
        dY(1+5*nodes:6*nodes)= (block.MethSpacing*QT(1+5*nodes:6*nodes) + Hin3 - Hout3)./tC(1+5*nodes:6*nodes);  %Fuel Reformer Channels
end
end%Ends function DynamicTemps