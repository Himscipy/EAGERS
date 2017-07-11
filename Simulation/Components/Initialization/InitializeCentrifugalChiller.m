function block = InitializeCentrifugalChiller(varargin)
%2 stage centrifugal chiller
% Three (3) Inlets: Cold water flow, cooling tower flow, Compressor power,
% Two (2) Outlets: Chilled water fow, coolint tower return flow
% Seven, (7) States: compressor RPM, 6 temperatures (evaporator: hot (Chilled water), solid, cold (refrigerant), and condensor: hot (refrigerant), solid, cold (CTwater))
global Tags
block = varargin{1};
Cp_H2O = 4.186;
if length(varargin)==1 % first initialization
    block.Area_ev = 2*block.Refrig_Tons*3.517/(block.hconv_ev*block.HeaterDeltaT);%200; %Surface area (m^2)
    block.Area_cd = 2*(block.Refrig_Tons*(1+1/block.EstimatedCOP))*3.517/(block.hconv_cd*block.HeaterDeltaT);%; %Surface area (m^2)
    
    Dir=strrep(which('InitializeCentrifugalChiller.m'),fullfile('Components','Initialization','InitializeCentrifugalChiller.m'),'CompressorMaps');
    load(fullfile(Dir,block.Map));
    f = fieldnames(map);
    for i = 1:1:length(f)
        block.(f{i}) = map.(f{i});
    end
    Eff = map.Efficiency;
    Eff(Eff<0)=nan;
    block.minEfficiency = min(min(Eff));
    block.Efficiency(block.Efficiency<0) = block.minEfficiency;
    
    block.Scale = [1e4*pi()/30 [4 3 2 30 29 28]]; %RPM & 6 temperatures
    block.IC = ones(length(block.Scale),1);
    %%
    block.InletPorts = {'Power','CWflow','CTflow'};
    block.CWflow.IC.T = 11+273.15;
    block.CWflow.IC.H2O = block.Refrig_Tons*3.517/(Cp_H2O*7)/18;
    block.CTflow.IC.T = 22+273.15;
    block.CTflow.IC.H2O = block.Refrig_Tons*(1+1/block.EstimatedCOP)*3.517/(Cp_H2O*5)/18;
    block.Power.IC = block.Refrig_Tons/block.EstimatedCOP;
    
    block.OutletPorts = {'ChilledWater';'CoolingTowerReturn';'ColdWaterTemp';'CTreturnTemp';};
    block.ChilledWater.IC.T = block.Scale(2)+273.15;
    block.ChilledWater.IC.H2O = block.Refrig_Tons*3.517/(Cp_H2O*(block.CWflow.IC.T - block.ChilledWater.IC.T));
    block.CoolingTowerReturn.IC.T = block.Scale(7)+273.15;
    block.CoolingTowerReturn.IC.H2O = block.Refrig_Tons*(1 + 1/block.EstimatedCOP)*3.517/(Cp_H2O*(block.CoolingTowerReturn.IC.T - block.CTflow.IC.T));
    block.ColdWaterTemp.IC = block.Scale(2);
    block.CTreturnTemp.IC = block.Scale(7);
    
    block.P_Difference = {};
    
end
if length(varargin)==2 %% Have inlets connected, re-initialize   
    Inlet = varargin{2};
    %need to solve problem, given power find temperature/pressure/speed states
    error = 1;
    Y = block.Scale;
    while abs(error)>.001
        %Start with guess of Y(4): T1. Find flow and PR for given Y(4)
        %find Qev & Y(2)   % in SS the evaporator HX temp shouldn't be changing
        error2 = 1;
        while abs(error2)>1e-3
            dT = Y(2)-Y(4);
            Qev = block.hconv_ev*block.Area_ev/2*dT; %heat transfer from water being cooled to refrigerant
            Q_C = (Inlet.CWflow.T-273.15 - Y(2))*(Inlet.CWflow.H2O*18*Cp_H2O);
            error2 = (Q_C/Qev-1); %error in heat transfer
            dT = dT*(1+.5*error2);
            Y(2) = Y(4) + dT;
        end
        Qcd = Qev + Inlet.Power*block.MotorEfficiency;
        %find Y(7) & Y(5) from condensor energy balance
        Y(7) = (Inlet.CTflow.T-273.15 + Qcd/(Inlet.CTflow.H2O*18*Cp_H2O));
        Y(5) = Qcd*2/(block.hconv_cd*block.Area_cd) + Y(7); %heat transfer from refrigerant to  cooling tower water loop
        
        Pev = refrigprop(block.Refrigerant,Y(4)+273.15,'T','P');% evaporator pressure
        Pcd = refrigprop(block.Refrigerant,Y(5)+273.15,'T','P');% condensor pressure
        PR = (Pcd/Pev);
        Pint = PR^.5*Pev; % multiply by .14504 and subtract 14.7 to get psig
        %Economizer
        H1 = refrigprop(block.Refrigerant,Y(4)+273.15,'T','Hsat_vap');
        H5 = refrigprop(block.Refrigerant,Y(5)+273.15,'T','Hsat_liq');
        H7 = refrigprop(block.Refrigerant,Pint,'P','Hsat_vap');
        H8 = refrigprop(block.Refrigerant,Pint,'P','Hsat_liq');
        X = (H5-H8)/(H7-H8); %quality from economizer?
        
        Flow = Qev/((1-X)*(H1-H8));
        
        % Find work associated with this Flow & PR
        error4 = 1;
        while abs(error4)>1e-3
            RPM = Y(1)*30/pi();
            Pnorm = (PR-1)/(block.compdesP-1) + 1; %feed into calculation of beta
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
                Beta = interp1(PRVec,block.Beta,Pnorm,'spline');
                Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
                Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
            end
            error4 = (Flow - Nflow*block.FlowDesign)/Flow;
            Y(1) = Y(1)*(1 + .33*error4);
        end
        %%Compressor
        % Work = Sat vapor @ Pev to mix in economizer at Pint, then all to Pcd
        S1 = refrigprop(block.Refrigerant,Pev,'P','Ssat_vap');
        H2s = refrigprop(block.Refrigerant,Pint,'P',S1,'S','H');
        H2 = (H2s - H1)/Eff + H1;

        H3 = H7*X + (1-X)*H2;
        S3 = refrigprop(block.Refrigerant,Pint,'P',H3,'H','S');
        H4s = refrigprop(block.Refrigerant,Pcd,'P',S3,'S','H');
        H4 = (H4s-H3)/Eff + H3;
        Work = Flow*((1-X)*(H2-H1) + (H4-H3));
                
        % Change Y(4) to match input work.
        error = (Inlet.Power*block.MotorEfficiency - Work)/Work;
        Y(4) = Y(4) - 5*error;
    end
    Y(3) = (Y(4) + Y(2))/2;
    Y(6) = (Y(7) + Y(5))/2;
    block.Scale = Y;
    Out.ChilledWater.T = Y(2) + 273.15; %(Inlet.CWflow.T - Qev/(Inlet.CWflow.H2O*18*Cp_H2O));
    Out.CoolingTowerReturn.T = Y(7) + 273.15; %(Inlet.CTflow.T + Qcd/(Inlet.CTflow.H2O*18*Cp_H2O));
    Out.ChilledWater.H2O = Inlet.CWflow.H2O;
    Out.CoolingTowerReturn.H2O = Inlet.CTflow.H2O;
    
    block.ChilledWater.IC = Out.ChilledWater;
    block.CoolingTowerReturn.IC = Out.CoolingTowerReturn;
    block.ColdWaterTemp.IC = Out.ChilledWater.T - 273.15;
    block.CTreturnTemp.IC = Out.CoolingTowerReturn.T - 273.15;
    Tags.(block.name).RPM = RPM;
    Tags.(block.name).EvaporatorPressure = Pev;
    Tags.(block.name).CondensorPressure = Pcd;
    Tags.(block.name).Cooling = Qev;
    Tags.(block.name).ColdWaterTemp = block.ColdWaterTemp.IC;
end