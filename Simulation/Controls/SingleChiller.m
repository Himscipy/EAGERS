function Out = SingleChiller(t,Y, Inlet,block,string1)
% Controller for single chiller plant
% Seven (7) Inlets: Cooling Power Demand, Chiller Tinlet, Chiller Toutlet, Chiller flow rate, CT Tinlet, CT Toutlet, CT flow rate
% Four (4) Outlets: Chiller power, water pump power, CT fan power, CT pump power, 
% Four (4) States: Chiller power, CT fan power, CW pump power, CT pump power
global Tags   
Cp_H2O = 4.186;
Chiller1Error = (Inlet.CW_Toutlet - block.Target(1))/(Inlet.CW_Tinlet - block.Target(1));% error in chiller exit temp, converted to a % error in electric power input.  flow converted from gallons per minute to kg/s * Cp
CTfanError = (Inlet.CT_Toutlet - block.Target(2))/3;% percent error in fan flow: assumes doubling/halving air flow would increase/decrease cooling 4C in cooling tower loop
idealCWflow = Inlet.CoolingPower/(Cp_H2O*(Inlet.CW_Tinlet - block.Target(1))); %mass flow that would result in correct cooling
CWpumpError = (idealCWflow - Inlet.CW_flow/15.83)/idealCWflow; %(Y(3)*(1+block.PropGain(3))/(1+block.PropGain(3)/idealCWflow))/idealCWflow;
idealCTflow = (Inlet.CoolingPower+block.Scale(1))/(Cp_H2O*(Inlet.CT_Tinlet - block.Target(2))); %mass flow that would result in correct cooling
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