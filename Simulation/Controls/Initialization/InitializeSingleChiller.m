function block = InitializeSingleChiller(varargin)
% Controller for single chiller plant
% Seven (7) Inlets: Cooling Power Demand, Chiller Tinlet, Chiller Toutlet, Chiller flow rate, CT Tinlet, CT Toutlet, CT flow rate
% Four (4) Outlets: Chiller power, water pump power, CT fan power, CT pump power, 
% Four (4) States: Chiller power, CT fan power, CW pump power, CT pump power
global Tags
block = varargin{1};
Cp_H2O = 4.186;
if length(varargin)==1 % first initialization
    block.description = {'Chiller power'; 'CT fan power';'CW pump power';'CT pump power';};
    
    block.InletPorts = {'CoolingPower','CW_Tinlet','CW_Toutlet','CW_flow','CT_Tinlet','CT_Toutlet','CT_flow'};
    block.CoolingPower.IC = block.NominalCoolingCapacity;
    block.CW_Tinlet.IC = 11;
    block.CW_Toutlet.IC = block.Target(1);
    block.CW_flow.IC = block.NominalCoolingCapacity/(Cp_H2O*(block.CW_Tinlet.IC - block.CW_Toutlet.IC))*18;
    block.CT_Tinlet.IC = 28;
    block.CT_Toutlet.IC = block.Target(2);
    block.CT_flow.IC = block.NominalCoolingCapacity/(Cp_H2O*(block.CT_Tinlet.IC - block.CT_Toutlet.IC))*18;
    
    block.OutletPorts = {'Chiller1';'CoolingTowerFan';'ColdWaterPump';'CoolingTowerPump';};
    block.Chiller1.IC = block.NominalCoolingCapacity/(block.EstimatedCOP+1);
    
%     AirFlow = (block.NominalCoolingCapacity+block.Chiller1.IC)/(1.006*6); %mass flow rate (kg/s of air  in cooling tower, assuming Cp = 1.006 kJ/kg and 6C temperature drop
    block.CoolingTowerFan.IC = 0.03*block.Chiller1.IC;% 3% of net power %AirFlow*9.81*0.3/(0.75*1e3);% fan power in kW, assuming 0.03m head loss(1 inch) , and 90% efficiency
    
%     WaterFlow = block.NominalCoolingCapacity/(Cp_H2O*7); %mass flow of cooling water (kg/s) assuming 7C of temperature drop
    block.ColdWaterPump.IC =  0.05*block.Chiller1.IC;% 5% of net power % WaterFlow*9.81*20/(0.75*1e3);% pump power in kW, assuming 3m head loss(5 psi), and 90% efficiency
    
%     CT_Flow = block.NominalCoolingCapacity/(Cp_H2O*6); %mass flow of cooling water (kg/s) assuming 6C of temperature drop
    block.CoolingTowerPump.IC =  0.04*block.Chiller1.IC;% 4% of net power %CT_Flow*9.81*10/(0.75*1e3);% pump power in kW, assuming 1m head loss(1.5 psi), and 90% efficiency
    
    block.Scale = [block.Chiller1.IC;block.CoolingTowerFan.IC;block.ColdWaterPump.IC;block.CoolingTowerPump.IC;]; 
    block.IC = ones(length(block.Scale),1);
    Tags.(block.name).COP =  block.NominalCoolingCapacity/sum(block.Scale);
    block.P_Difference = {};
end 
if length(varargin)==2 %% Have inlets connected, re-initialize  
    Inlet = varargin{2};
    idealCWflow = Inlet.CoolingPower/(Cp_H2O*(Inlet.CW_Tinlet - block.Target(1))); %mass flow that would result in correct cooling
    idealCTflow = (Inlet.CoolingPower+block.Scale(1))/(Cp_H2O*(Inlet.CT_Tinlet - block.Target(2))); %mass flow that would result in correct cooling
    
    error =[(Inlet.CW_Toutlet - block.Target(1))*Inlet.CW_flow/15.83*Cp_H2O/Tags.(block.name).COP;...;%flow converted from gallons per minute to kg/s * Cp
            block.Scale(2)*(Inlet.CT_Toutlet - block.Target(2))/3;... % assumes doubling/halving air flow would increase/decrease cooling 3C in cooling tower loop
            block.Scale(3)*(idealCWflow - Inlet.CW_flow/15.83)/idealCWflow;... %convert GPM to kg/s 
            block.Scale(4)*(idealCTflow - Inlet.CT_flow/15.83)/idealCTflow;]; %convert GPM to kg/s 
    errorNorm = error./block.Scale;
    block.Scale = block.Scale + error;
    
    block.Chiller1.IC = block.Scale(1);
    block.CoolingTowerFan.IC = block.Scale(2);
    block.ColdWaterPump.IC = block.Scale(3);
    block.CoolingTowerPump.IC = block.Scale(4);
    
    block.InitializeError = max(abs(errorNorm));
    
    Tags.(block.name).Chilling = (Inlet.CW_Tinlet - block.Target(1))*Inlet.CW_flow/15.83*Cp_H2O;
    Tags.(block.name).COP =  Tags.(block.name).Chilling/sum(block.Scale);
    Tags.(block.name).Power = sum(block.Scale); %power in kW
end