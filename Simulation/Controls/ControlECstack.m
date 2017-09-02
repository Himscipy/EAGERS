function Out = ControlECstack(varargin)
% Controls for electrolyzer stack only, control air flow to anode, inlet temperature and steam flow rate
% Four (4) inlets: Two Targets (temperature set point and power), Measured Temperature, Voltage
% Seven (7) outlets: Two measurements of targets, OxidentTemp, oxidant flow rate, steam temperature, steam flow rate, Current
% One  (1) state: If there is no oxidant flow, then control current to maintain temperature
% Like all controllers, if there are n targets, the first n inlets set 
% those targets, while the first n outlets correspond to the controllers 
% ability to hit that target
% The model can thus be linearized around any of these targets
global Tags
F=96485.339; % %Faraday's constant in Coulomb/mole
if length(varargin)==1 %first initialization
    block = varargin{1};
    block.PIdescription = {'Oxidant Flow Rate';'Oxidant Temperature';}; % THIS IS MY CHANGE
    block.TargetDescription = {'Operating Temperature';'Net Power'};
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = ComponentProperty(block.Target{j});%find the value of the desired properties in block and assing it to Target(j)
    end
    block.Target  = Target;%ovverwrite block.Target with the numerical values
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};%don't know why end+1
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);%create new field in block for the target
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);%and the measured value
    end
    block.InletPorts(end+1:end+2) = {'Temperature','Voltage'};%add two inlet ports for Temperature and Voltage
    block.Temperature.IC = ComponentProperty(strcat(block.connections{3},'.IC'));%looking for the value of EC1.MeasureTpen.IC
    block.Voltage.IC = ComponentProperty(strcat(block.connections{4},'.IC'));%looking for the value of EC1.MeasureVoltage.IC
    
    block.deltaTStack = ComponentProperty(block.deltaTStack);
    block.Cells = ComponentProperty(block.Cells);
    block.Steam = ComponentProperty(block.Steam);
    block.Utilization = ComponentProperty(block.Utilization);
    block.SteamTemperature = ComponentProperty(block.SteamTemperature);
    
    OxFlow = NetFlow(ComponentProperty(block.OxidantFlow)); 
    OxTemp = 1023;
    Current = block.Target2.IC*1000/(block.Voltage.IC*block.Cells);
    SteamFlow = (block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    [h,~] = enthalpy(block.Target1.IC,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    block.Vbalance = 1./(2*F)*h_rxn3; %voltage that balances heat

    block.OutletPorts = {'OxidantTemp','OxidantFlow','SteamTemp','SteamFlow','Current'};
    block.OxidantTemp.IC = block.Target1.IC;
    block = rmfield(block,'OxidantFlow'); % THIS IS MY CHANGE
    block.OxidantFlow.IC = OxFlow; 
    block.SteamTemp.IC = block.SteamTemperature;
    block.SteamFlow.IC = SteamFlow;
    block.Current.IC = Current;
    
    block.P_Difference = {};

    block.IC = [1;1]; % inital condition
        
    if OxFlow>0
        block.HasFlow = true;
        block.Scale = [OxFlow;OxTemp];
        block.UpperBound = [inf,inf];
        block.LowerBound = [0,0];
    else 
        block.HasFlow = false;
        block.Scale = [Current;];
        block.UpperBound = inf;
        block.LowerBound = -inf;
    end
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    block = varargin{1};
    Inlet = varargin{2};
    PowerError = (Inlet.Target2-abs(block.Current.IC)*Inlet.Voltage*block.Cells/1000)/Inlet.Target2;
    block.Current.IC = block.Current.IC*(1 + PowerError);
    block.SteamFlow.IC = (block.Cells*abs(block.Current.IC)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    block.SteamTemp.IC = block.SteamTemp.IC; %currently not manipulated during initialization
    
    averageT = mean(Inlet.Temperature);
    
    block.Measured1.IC = averageT;
    block.Measured2.IC = block.Current.IC*Inlet.Voltage*block.Cells/1000; %Power in kW
    
%     if block.HasFlow
%         Q_cathode = block.OxidantFlow.IC*40*block.deltaTStack;
%         
%         if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(block.Current.IC)/1000) - Q_cathode)>0
%             block.OxidantTemp.IC = Inlet.Target1-50; %cooling stack
%             TavgError = (averageT -Inlet.Target1)/block.deltaTStack; % too hot = increase flow
%         else
%             block.OxidantTemp.IC = Inlet.Target1+50;%heating stack
%             TavgError = (Inlet.Target1-averageT)/block.deltaTStack; %too hot = reduce flow
%         end
%         block.InitializeError = abs(TavgError); 
%         %newTemp=block.OxidantTemp.IC*(1+.5*TavgError);
%         %block.OxidantTemp.IC = newTemp;
%         %block.OxidantTemp.IC = 923;
%         newFlow = block.OxidantFlow.IC*(1+.5*TavgError);
%         block.OxidantFlow.IC = newFlow;
%         block.Scale = [block.OxidantFlow.IC,1023];
%         block.IC = [1-TavgError*block.PropGain(1),1-TavgError*block.PropGain(2)]; %Oxidant flow rate
%     else
%         block.InitializeError = 0;
%         block.Scale = [block.Current.IC];
%         block.IC = [1-PowerError*block.PropGain(2)]; % inital condition with no anode flow in: net current
%     end

    if block.HasFlow
        if ~isfield (Tags.EC1,'StackdeltaT')
            Q_cathode = block.OxidantFlow.IC*40*block.deltaTStack;
        else
             Q_cathode = block.OxidantFlow.IC*SpecHeat(Tags.AirSource.Outlet)* Tags.EC1.StackdeltaT;
        end
        
        if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(block.Current.IC)/1000))>0
            if ~isfield (Tags.EC1,'StackdeltaT')
            block.OxidantTemp.IC = Inlet.Target1-80; % Change this number (60) if crashes near thermoneutral
            newTemp=block.OxidantTemp.IC;
            newFlow=block.OxidantFlow.IC;
            TavgError = (averageT -Inlet.Target1)/block.deltaTStack;
            else
            
            TavgError = (averageT -Inlet.Target1)/block.deltaTStack; 
            if Q_cathode>0
                newTemp=block.OxidantTemp.IC*(1-.005*TavgError);
                newFlow=block.OxidantFlow.IC*(1+.2*TavgError);
                
            else
                newTemp=block.OxidantTemp.IC*(1-.005*TavgError);
                newFlow=block.OxidantFlow.IC*(1-.2*TavgError);
                
            end
            end
        else
            if ~isfield (Tags.EC1,'StackdeltaT')
            block.OxidantTemp.IC = Inlet.Target1+80; % Change this number (60) if crashes near thermoneutral
            newTemp=block.OxidantTemp.IC;
            newFlow=block.OxidantFlow.IC;
            TavgError = (Inlet.Target1-averageT)/block.deltaTStack;
            else
            
            TavgError = (Inlet.Target1-averageT)/block.deltaTStack; 
            if Q_cathode>0
                newTemp=block.OxidantTemp.IC*(1+.005*TavgError);
                newFlow=block.OxidantFlow.IC*(1-.2*TavgError);
            else
                newTemp=block.OxidantTemp.IC*(1+.005*TavgError);
                newFlow=block.OxidantFlow.IC*(1+.2*TavgError);
            end
            end
        end
        block.InitializeError = abs(TavgError); 
        block.OxidantFlow.IC = newFlow;
        block.OxidantTemp.IC = newTemp;
        block.Scale = [block.OxidantFlow.IC,block.OxidantTemp.IC];
        if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(block.Current.IC)/1000))>0
            if Q_cathode>0
                block.IC = [1+TavgError*block.PropGain(1),1-TavgError*block.PropGain(2)]; %Oxidant flow rate
            else
                 block.IC = [1-TavgError*block.PropGain(1),1-TavgError*block.PropGain(2)];
            end
        else
            if Q_cathode>0
                 block.IC = [1-TavgError*block.PropGain(1),1+TavgError*block.PropGain(2)];
            else
                block.IC = [1+TavgError*block.PropGain(1),1+TavgError*block.PropGain(2)];
            end
        end
    else
        block.InitializeError = 0;
        block.Scale = [block.Current.IC];
        block.IC = [1-PowerError*block.PropGain(2)]; % inital condition with no anode flow in: net current
    end
    Out = block;
else %running the model
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    averageT = mean(Inlet.Temperature);
    % averageT = (Inlet.Hot - Inlet.Cold + block.dT_cath_PEN);

    if block.HasFlow
        Current = -Inlet.Target2*1000/(Inlet.Voltage*block.Cells);
    else
        Current = Y(1);
        VoltError = block.Vbalance - Inlet.Voltage;
    end
    SteamFlow = block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000;

    if block.HasFlow
        
         Q_cathode = Y(1)*SpecHeat(Tags.AirSource.Outlet)* Tags.EC1.StackdeltaT;   % THIS IS MY CHANGE
        %Q_cathode = Y(1)*40*block.deltaTStack;  
        
%         if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(Current)/1000) - Q_cathode)>0
%             TavgError = (averageT-Inlet.Target1)/block.deltaTStack; % too hot = increase flow
%             
%                 Out.OxidantTemp = min(1123,Y(2)-(TavgError*1*block.PropGain(2))*block.Scale(2));
%                 Out.OxidantTemp = max(923,Out.OxidantTemp);
%                 
%             
%         else
%             TavgError = (Inlet.Target1-averageT)/block.deltaTStack; %too hot = reduce flow
%             
%                 Out.OxidantTemp = min(1123,Y(2)+(TavgError*1*block.PropGain(2))*block.Scale(2));
%                 Out.OxidantTemp = max(923,Out.OxidantTemp);
%                 
%             
%         end
%         Out.OxidantFlow = max(0,Y(1)+(TavgError*block.PropGain(1))*block.Scale(1));
        
        

        if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(Current)/1000))>0  % THIS IS MY CHANGE
            
            TavgError = (averageT-Inlet.Target1)/block.deltaTStack; 
           
            if Q_cathode>0
                Out.OxidantTemp = min(1123,Y(2)-(TavgError*1*block.PropGain(2))*block.Scale(2));  % THIS IS MY CHANGE
                Out.OxidantTemp = max(923,Out.OxidantTemp);  % THIS IS MY CHANGE
                Out.OxidantFlow = max(0,Y(1)+(TavgError*block.PropGain(1))*block.Scale(1));  % THIS IS MY CHANGE
            else
                Out.OxidantTemp = min(1123,Y(2)-(TavgError*1*block.PropGain(2))*block.Scale(2));  % THIS IS MY CHANGE
                Out.OxidantTemp = max(923,Out.OxidantTemp);  % THIS IS MY CHANGE
                Out.OxidantFlow = max(0,Y(1)-(TavgError*block.PropGain(1))*block.Scale(1));  % THIS IS MY CHANGE
            end
                
            
        else
            
            TavgError = (Inlet.Target1-averageT)/block.deltaTStack; 
            
            if Q_cathode>0
                Out.OxidantTemp = min(1123,Y(2)+(TavgError*1*block.PropGain(2))*block.Scale(2));  % THIS IS MY CHANGE
                Out.OxidantTemp = max(923,Out.OxidantTemp);  % THIS IS MY CHANGE
                Out.OxidantFlow = max(0,Y(1)-(TavgError*block.PropGain(1))*block.Scale(1));  % THIS IS MY CHANGE

            else
                Out.OxidantTemp = min(1123,Y(2)+(TavgError*1*block.PropGain(2))*block.Scale(2));  % THIS IS MY CHANGE
                Out.OxidantTemp = max(923,Out.OxidantTemp);  % THIS IS MY CHANGE
                Out.OxidantFlow = max(0,Y(1)+(TavgError*block.PropGain(1))*block.Scale(1));  % THIS IS MY CHANGE
                
            end
                
            
        end
       
    else
        Out.OxidantTemp = Inlet.Target1;
        Out.OxidantFlow = 0;
    end

    if strcmp(string1,'Outlet')
        Out.Measured1.IC = averageT;
        Out.Measured2.IC = Current*Inlet.Voltage*block.Cells/1000; %Power in kW
        Out.SteamTemp = block.SteamTemp.IC;%currently no control of fuel inlet temp
        Out.SteamFlow = SteamFlow;
        Out.Current = Current;
        Tags.(block.name).OxidantTemp = Out.OxidantTemp;
        Tags.(block.name).OxidantFlow = Out.OxidantFlow;
        Tags.(block.name).Tsteam = block.SteamTemp.IC;%currently no control of steam inlet temp
        Tags.(block.name).SteamFlow = SteamFlow;
        Tags.(block.name).Current = Current;
    elseif strcmp(string1,'dY')
        dY = Y*0;
        if block.HasFlow
            if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(Current)/1000))>0
                if Q_cathode>0
               
                   dY(2) = -block.Gain(2)*TavgError;  % THIS IS MY CHANGE
                   dY(1) = block.Gain(1)*TavgError;   % THIS IS MY CHANGE
                else
                    dY(2) = -block.Gain(2)*TavgError; % THIS IS MY CHANGE
                   dY(1) = -block.Gain(1)*TavgError;  % THIS IS MY CHANGE
                    
                end
            else
                if Q_cathode>0
               
                   dY(2) = block.Gain(2)*TavgError;  % THIS IS MY CHANGE
                   dY(1) = -block.Gain(1)*TavgError;  % THIS IS MY CHANGE
                else
                    dY(2) = block.Gain(2)*TavgError;  % THIS IS MY CHANGE
                    dY(1) = block.Gain(1)*TavgError;  % THIS IS MY CHANGE
                    
                end
            end
            
        else
            dY(1) = block.Gain(1)*VoltError;
        end
        
        
        %         if block.HasFlow
%             if ((block.Cells*(Inlet.Voltage - block.Vbalance)*abs(Current)/1000) - Q_cathode)>0
%                
%                    dY(2) = -block.Gain(2)*TavgError;
%             else
%                 
%                    dY(2) = block.Gain(2)*TavgError;
%             end
%             if Y(1)>0 || TavgError>0 %anti-windup so flow does not go negative
%                 dY(1) = block.Gain(1)*TavgError;
%                 
%             end
%         else
%             dY(1) = block.Gain(1)*VoltError;
%         end
        Out = dY;
    end
end
end%Ends function ControlECstack
