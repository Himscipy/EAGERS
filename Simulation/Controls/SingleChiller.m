function Out = SingleChiller(varargin)
% Controller for single chiller plant with a cooling tower loop
% Nine (9) Inlets: Three targets (Cooling Power Demand, Chilled water temp, cooling tower return temp), Chiller Tinlet, Chiller Toutlet, Chiller flow rate, CT Tinlet, CT Toutlet, CT flow rate
% Seven (7) Outlets: Three associated with the targets, Chiller power, water pump power, CT fan power, CT pump power, 
% Four (4) States: Chiller power, CT fan power, CW pump power, CT pump power
global Tags   
Cp_H2O = 4.186;
if length(varargin)==1 % first initialization
    block = varargin{1};
    block.PIdescription = {'Chiller power'; 'CT fan power';'CW pump power';'CT pump power';};
    block.TargetDescription = {'Cooling Power Demand';'Chilled Water Temperature';'Cooling Tower Return Temperature';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(block.Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target = Target;
    
    block.InletPorts = {};
    block.OutletPorts = {};
    for i = 1:1:length(block.Target)
        block.InletPorts(end+1) = {strcat('Target',num2str(i))};
        block.OutletPorts(end+1) = {strcat('Measured',num2str(i))};
        block.(strcat('Target',num2str(i))).IC = block.Target(i);
        block.(strcat('Measured',num2str(i))).IC = block.Target(i);
    end
    block.InletPorts(end+1:end+6) = {'CW_Tinlet','CW_Toutlet','CW_flow','CT_Tinlet','CT_Toutlet','CT_flow'};
    block.CW_Tinlet.IC = 11;
    block.CW_Toutlet.IC = block.Target2.IC;
    block.CW_flow.IC = block.Target1.IC/(Cp_H2O*(block.CW_Tinlet.IC - block.CW_Toutlet.IC))*18;
    block.CT_Tinlet.IC = 28;
    block.CT_Toutlet.IC = block.Target3.IC;
    block.CT_flow.IC = block.Target1.IC/(Cp_H2O*(block.CT_Tinlet.IC - block.CT_Toutlet.IC))*18;
    block.OutletPorts(end+1:end+4) = {'Chiller1';'CoolingTowerFan';'ColdWaterPump';'CoolingTowerPump';};

    block.Chiller1.IC = block.Target1.IC/(block.EstimatedCOP+1);
    
%     AirFlow = (block.Target1.IC+block.Chiller1.IC)/(1.006*6); %mass flow rate (kg/s of air  in cooling tower, assuming Cp = 1.006 kJ/kg and 6C temperature drop
    block.CoolingTowerFan.IC = 0.03*block.Chiller1.IC;% 3% of net power %AirFlow*9.81*0.3/(0.75*1e3);% fan power in kW, assuming 0.03m head loss(1 inch) , and 90% efficiency
    
%     WaterFlow = block.Target1.IC/(Cp_H2O*7); %mass flow of cooling water (kg/s) assuming 7C of temperature drop
    block.ColdWaterPump.IC =  0.05*block.Chiller1.IC;% 5% of net power % WaterFlow*9.81*20/(0.75*1e3);% pump power in kW, assuming 3m head loss(5 psi), and 90% efficiency
    
%     CT_Flow = block.Target1.IC/(Cp_H2O*6); %mass flow of cooling water (kg/s) assuming 6C of temperature drop
    block.CoolingTowerPump.IC =  0.04*block.Chiller1.IC;% 4% of net power %CT_Flow*9.81*10/(0.75*1e3);% pump power in kW, assuming 1m head loss(1.5 psi), and 90% efficiency
    
    block.Scale = [block.Chiller1.IC;block.CoolingTowerFan.IC;block.ColdWaterPump.IC;block.CoolingTowerPump.IC;]; 
    block.IC = ones(length(block.Scale),1);
    Tags.(block.name).COP =  block.Target1.IC/sum(block.Scale);
    block.P_Difference = {};
    Out = block;
elseif length(varargin)==2 %% Have inlets connected, re-initialize  
    block = varargin{1};
    Inlet = varargin{2};
    idealCWflow = Inlet.Target1/(Cp_H2O*(Inlet.CW_Tinlet - Inlet.Target2)); %mass flow that would result in correct cooling
    idealCTflow = (Inlet.Target1+block.Scale(1))/(Cp_H2O*(Inlet.CT_Tinlet - Inlet.Target3)); %mass flow that would result in correct cooling
    
    error =[(Inlet.CW_Toutlet - Inlet.Target2)*Inlet.CW_flow/15.83*Cp_H2O/Tags.(block.name).COP;...;%flow converted from gallons per minute to kg/s * Cp
            block.Scale(2)*(Inlet.CT_Toutlet - Inlet.Target3)/3;... % assumes doubling/halving air flow would increase/decrease cooling 3C in cooling tower loop
            block.Scale(3)*(idealCWflow - Inlet.CW_flow/15.83)/idealCWflow;... %convert GPM to kg/s 
            block.Scale(4)*(idealCTflow - Inlet.CT_flow/15.83)/idealCTflow;]; %convert GPM to kg/s 
    errorNorm = error./block.Scale;
    block.Scale = block.Scale + error;
    
    block.Chiller1.IC = block.Scale(1);
    block.CoolingTowerFan.IC = block.Scale(2);
    block.ColdWaterPump.IC = block.Scale(3);
    block.CoolingTowerPump.IC = block.Scale(4);
    
    block.InitializeError = max(abs(errorNorm));
    
    Tags.(block.name).Chilling = (Inlet.CW_Tinlet - Inlet.Target2)*Inlet.CW_flow/15.83*Cp_H2O;
    Tags.(block.name).COP =  Tags.(block.name).Chilling/sum(block.Scale);
    Tags.(block.name).Power = sum(block.Scale); %power in kW
    Out = block;
else
    t = varargin{1};
    Y = varargin{2};
    Inlet = varargin{3};
    block = varargin{4};
    string1 = varargin{5};
    Chiller1Error = (Inlet.CW_Toutlet - Inlet.Target2)/(Inlet.CW_Tinlet - Inlet.Target2);% error in chiller exit temp, converted to a % error in electric power input.  flow converted from gallons per minute to kg/s * Cp
    CTfanError = (Inlet.CT_Toutlet - Inlet.Target3)/3;% percent error in fan flow: assumes doubling/halving air flow would increase/decrease cooling 4C in cooling tower loop
    idealCWflow = Inlet.Target1/(Cp_H2O*(Inlet.CW_Tinlet - Inlet.Target2)); %mass flow that would result in correct cooling
    CWpumpError = (idealCWflow - Inlet.CW_flow/15.83)/idealCWflow; %(Y(3)*(1+block.PropGain(3))/(1+block.PropGain(3)/idealCWflow))/idealCWflow;
    idealCTflow = (Inlet.Target1+block.Scale(1))/(Cp_H2O*(Inlet.CT_Tinlet - Inlet.Target3)); %mass flow that would result in correct cooling
    CTpumpError = (idealCTflow - Inlet.CT_flow/15.83)/idealCTflow; %(Y(4)*(1+block.PropGain(4))/(1+block.PropGain(4)/idealCTflow))/idealCTflow;
    if strcmp(string1,'Outlet')
        Out.Chiller1 = Y(1)*(1 + Chiller1Error*block.PropGain(1));
        Out.CoolingTowerFan = Y(2)*(1 + CTfanError*block.PropGain(2));
        Out.ColdWaterPump = Y(3)*(1 + CWpumpError*block.PropGain(3));
        Out.CoolingTowerPump = Y(4)*(1 + CTpumpError*block.PropGain(4));
        Tags.(block.name).Chilling = (Inlet.CW_Tinlet - Inlet.CW_Toutlet)*Inlet.CW_flow/15.83*Cp_H2O;
        Tags.(block.name).Power = (Out.Chiller1 + Out.CoolingTowerFan + Out.ColdWaterPump + Out.CoolingTowerPump); %power in kW
        Tags.(block.name).COP =  Tags.(block.name).Chilling/Tags.(block.name).Power;
    elseif strcmp(string1,'dY')
        dY(1) = Chiller1Error*block.Gain(1);
        dY(2) = CTfanError*block.Gain(2);
        dY(3) = CWpumpError*block.Gain(3);
        dY(4) = CTpumpError*block.Gain(4);
        Out = dY;
    end
end
end%Ends function SingleChiller