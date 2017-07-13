function Out = CentrifugalChiller(t,Y, Inlet,block,string1)
%2 stage centrigugal chiller
% Inlets: Cold water flow, cooling tower flow, Compressor power,
% Outlets: Chilled water fow, coolint tower return flow
% States: compressor RPM, 6 temperatures (evaporator: hot, solid, cold, and condensor: hot, solid, cold)
global Tags    
Cp_H2O = 4.186;
%%Evaporator   
Qhot_ev = block.hconv_ev*block.Area_ev*(Y(2)-Y(3)); %heat transfer from water being cooled to evaporator solid
Qcold_ev = block.hconv_ev*block.Area_ev*(Y(3)-Y(4)); %heat transfer from evaporator solid to refrigerant
%%Condensor
Qhot_cd = block.hconv_cd*block.Area_cd*(Y(5)-Y(6)); %heat transfer from refrigerant to condensor solid
Qcold_cd = block.hconv_cd*block.Area_cd*(Y(6)-Y(7)); %heat transfer from condensor solid to cooling tower water loop
Out.ChilledWater.T = (Inlet.CWflow.T - Qhot_ev/(Inlet.CWflow.H2O*18*Cp_H2O));
Out.CoolingTowerReturn.T = (Inlet.CTflow.T + Qcold_cd/(Inlet.CTflow.H2O*18*Cp_H2O));
Out.ColdWaterTemp = Out.ChilledWater.T - 273.15;
Out.CTreturnTemp = Out.CoolingTowerReturn.T - 273.15;
RPM = Y(1)*30/pi();
Pev = refrigprop(block.Refrigerant,Y(4)+273.15,'T','P');% evaporator pressure
Pcd = refrigprop(block.Refrigerant,Y(5)+273.15,'T','P');% condensor pressure,
Tags.(block.name).RPM = RPM;
Tags.(block.name).EvaporatorPressure = Pev;
Tags.(block.name).CondensorPressure = Pcd;
Tags.(block.name).Cooling = Qhot_ev;
if strcmp(string1,'Outlet')
    Out.ChilledWater.H2O = Inlet.CWflow.H2O;
    Out.CoolingTowerReturn.H2O = Inlet.CTflow.H2O;
elseif strcmp(string1,'dY')
    
    %%CompressorMap
    PR = (Pcd/Pev);
    Pint = PR^.5*Pev; % multiply by .14504 and subtract 14.7 to get psig
    Pnorm = (PR-1)/(block.compdesP-1) + 1; %feed into calculation of beta

    %compressor map
    nRPM =RPM/block.RPMdesign;%normalized RPM
    i2 = find(block.RPM >=nRPM,1,'first');
    if isempty(i2)
        i2 = length(block.RPM);
    elseif i2==1
        i2 =2;
    end
    i1 = i2-1;
    s =(nRPM - block.RPM(i1))/(block.RPM(i2) - block.RPM(i1));
    PRVec = block.PressRatio(i1,:)*(1-s) + block.PressRatio(i2,:)*s;
    if Pnorm>PRVec(end)
    %     disp('Stall Reached')
        Beta = Pnorm/PRVec(end);
        Eff = block.minEfficiency;
        flow_RPM = block.NflowGMap(i1,end)*(1-s) + block.NflowGMap(i2,end)*s;
        Nflow = 1/Beta*flow_RPM;
    else
        if any(isnan(PRVec)) || isnan(Pnorm) || isinf(Pnorm)
            disp('WTF')
        end
        Beta = interp1(PRVec,block.Beta,Pnorm,'spline');
        Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
        Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
    end
    Flow = Nflow*block.FlowDesign;
    %Economizer
    H1 = refrigprop(block.Refrigerant,Y(4)+273.15,'T','Hsat_vap');
    H5 = refrigprop(block.Refrigerant,Y(5)+273.15,'T','Hsat_liq');
    H7 = refrigprop(block.Refrigerant,Pint,'P','Hsat_vap');
    H8 = refrigprop(block.Refrigerant,Pint,'P','Hsat_liq');
    X = (H5-H8)/(H7-H8); %quality

    %%Compressor
    % Work = liquid compression and vapor compression
    S1 = refrigprop(block.Refrigerant,Pev,'P','Ssat_vap');
    H2s = refrigprop(block.Refrigerant,Pint,'P',S1,'S','H');
    H2 = (H2s - H1)/Eff + H1;

    H3 = H7*X + (1-X)*H2;
    S3 = refrigprop(block.Refrigerant,Pint,'P',H3,'H','S');
    H4s = refrigprop(block.Refrigerant,Pcd,'P',S3,'S','H');
    H4 = (H4s-H3)/Eff + H3;

    Work = Flow*((1-X)*(H2-H1) + (H4-H3));
    Cp_R123 = refrigprop(block.Refrigerant,Y(4)+273.15,'T','Cp');
    
    %change in states
    dY = zeros(length(block.IC),1);
    MoI = block.ShaftDensity*block.ShaftLength*pi()*block.ShaftRadius^4;
    if RPM<=1.2*block.RPMdesign && RPM>0.5*block.RPMdesign
        dY(1) = (Inlet.Power*block.MotorEfficiency - Work)*1000/(MoI*Y(1));
    end
    dY(2) = Out.ChilledWater.T - (Y(2)+273.15);
    dY(3) = (Qhot_ev-Qcold_ev)/(block.Mass_ev*block.CP_ev);
    dY(4) = (Qcold_ev + H8*(Qcold_ev/(H1-H8)) - (1-X)*Flow*H1)/(block.Volume_ev*22.4*Cp_R123);
    dY(5) = (Flow*H4 - Qhot_cd - H5*(Qhot_cd/(H4-H5)))/(block.Volume_cd*22.4*Cp_R123);
    dY(6) = (Qhot_cd-Qcold_cd)/(block.Mass_cd*block.CP_cd);
    dY(7) = Out.CoolingTowerReturn.T - (Y(7)+273.15);
    Out = dY;
end