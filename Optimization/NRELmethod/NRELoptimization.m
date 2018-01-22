function GenDisp = NRELoptimization(IC,Forecast,scaleCost)
global Plant
% Specifically set up for a plant with 1 electric utility, 1 FC, 1 chiller 1 battery in that order
dt = Plant.optimoptions.Resolution*3600;% time step for discrete model, seconds
nS = Plant.optimoptions.Horizon/Plant.optimoptions.Resolution;% # of time steps in optimization
N_schedule = Plant.optimoptions.Horizon; % assumes hourly scheduling
%%need to replace these with inputs from GUI stored in Plant structure
% cas = 0.1*ones(nS,1); % ancillary service price
% cas(1:nS/2) = 0.001;
% cas(nS*3/4:end) = 0.001;
% 
% Forecast.Demand.H = 10*ones(nS,1); % incidental heating load in kW
% Forecast.Demand.H(1:nS/4) = 5;
% Forecast.Demand.H(nS*3/4:end) = 5;
Forecast.Demand.H = Forecast.Demand.H*2;

% T_high = 22.78; % temperature bounds
% T_low = 21.67;
T_high = 25; % temperature bounds
T_low = 20;
kp = -10; % controller gains
ki = -0.01;
alpha = 0.4; % outdoor air ratio
Cp = 1e3; % specific heat of air, J/(kg Kelvin)
Cr = 5e7; % room thermal capacitance, J/Kelvin
Cw = 5e7; % wall thermal capacitance, J/Kelvin
R = 1e-3; % wall thermal resistance, Kelvin/Watt
Ts = 12.78; % supply air temperature (celcius)
tau_ch = 200; % time constant of the chiller, second
tau_fc = 120*60; % time constant of the FC, second
alpha_f = 500; % fan power coefficient
T_approx = 22.22; % approximation for room temperature (celcius)
Tref_high = 38; % temperature setpoint bounds (celcius)
Tref_low = -20;%(celcius)
% kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service

eta_fce = Plant.Generator(3).Output.Electricity(end); % electrical efficiency of the FC
Pfcf_max = Plant.Generator(3).Size*1000; % FC capacity (W)
COP = Plant.Generator(4).Output.Cooling(end); %chiller efficiency

DischCurrent = Plant.Generator(5).VariableStruct.PeakDisch.*Plant.Generator(5).Size/Plant.Generator(5).VariableStruct.Voltage*1000;
DischResistScaled = (100/DischCurrent)*Plant.Generator(5).VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
DischVoltLoss = DischCurrent*DischResistScaled; %keep in mind when calculating loss as function of discharge current
ChargeCurrent = Plant.Generator(5).VariableStruct.PeakCharge*Plant.Generator(5).Size/Plant.Generator(5).VariableStruct.Voltage*1000;
ChargeResistScaled = (100/ChargeCurrent)*Plant.Generator(5).VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
ChargeVoltLoss = ChargeCurrent*ChargeResistScaled;
eta_c = Plant.Generator(5).VariableStruct.Voltage/(Plant.Generator(5).VariableStruct.Voltage+ChargeVoltLoss); %charging efficiency
eta_d = (Plant.Generator(5).VariableStruct.Voltage-DischVoltLoss)/Plant.Generator(5).VariableStruct.Voltage; %discharging efficiency
Pbat_max = (DischCurrent*Plant.Generator(5).VariableStruct.Voltage*eta_d); % battery ramping capacity (J/s)
Ebat_max = Plant.Generator(5).Size*(Plant.Generator(5).VariableStruct.MaxDOD/100)*1000*3600; % usable battery capacity (J)

%% Determine the IC

% 6 states: [T,Tw,mpi,Pch,Pfcf,Ebat]
% [room temp, wall temp, PI controller state, chiller power, FC fuel consumption, battery SOC]
% 4 schedule variables: [T_ref,Pfcf_ref,Pbatin,Pbatout]
% [room temp setpoint, FC fuel consumption setpoint, battery discharging, battery charging]

T0 = 22.22;
Tw0 = 26.667;
mpi0 = 0;
Pfcf0 = IC(3)*1000/eta_fce;%fuel cell fuel input in W
Pch0 = IC(4)*1000;%chiller cooling power in W
Ebat0 = IC(5)*1000*3600;%battery state of charge in J

Ti = T0;
Twi = Tw0;
mpii = mpi0;
Pchi = Pch0;
% Pfcfi = Pfcf0;
% Ebati = Ebat0;

T_refi = T0;
% Pfcf_refi = Pfcf0;
% Pbat_discharge_i = 0;
% Pbat_charge_i = 0;

for i = 1:nS
    erri = T_refi - Ti;
    mi = kp*erri + ki*mpii;
    Tmixi = alpha*Forecast.Temperature(1) + (1-alpha)*T_approx;
    Pcci = mi*Cp*(Tmixi-Ts);
    
    T_new = Ti + dt*1/Cr*(1/R*(Twi-Ti) + mi*Cp*(Ts-T_approx) + Forecast.Demand.H(1)*1000);
    Tw_new = Twi + dt*1/Cw*(1/R*(Ti-Twi) + 1/R*(Forecast.Temperature(1)-Twi));
    mpi_new = mpii + dt*erri;
    Pch_new = Pchi + dt*1/tau_ch*(-Pchi + Pcci);
%     Pfcf_new = Pfcfi + dt*1/tau_fc*(-Pfcfi + Pfcf_refi);
%     Ebat_new = Ebati - dt*1/eta_d*Pbat_discharge_i + dt*eta_c*Pbat_charge_i;
    
    Ti = T_new;
    Twi = Tw_new;
    mpii = mpi_new;
    Pchi = Pch_new;
%     Pfcfi = Pfcf_new;
%     Ebati = Ebat_new;  
end

Tw0 = Twi;
mpi0 = mpii;



%% Opt

% cvx_precision low
cvx_begin
% cvx_solver_settings( 'dumpfile', 'test' ) 
cvx_solver_settings( 'NumericFocus', '3' ) 

    variables x(6*nS) u(4*N_schedule) %R_as(N_schedule)
    
    for i = 1:nS
        
        hour = ceil(i*Plant.optimoptions.Resolution);
   
        Ti = x((i-1)*6+1);
        Twi = x((i-1)*6+2);
        mpii = x((i-1)*6+3);
        Pchi = x((i-1)*6+4);
        Pfcfi = x((i-1)*6+5);
        Ebati = x((i-1)*6+6);

        T_refi = u((hour-1)*4+1);
        Pfcf_refi = u((hour-1)*4+2);
        Pbat_discharge_i = u((hour-1)*4+3);
        Pbat_charge_i = u((hour-1)*4+4);
        
%         R_asi = R_as(hour);

        erri = T_refi - Ti;
        mi = kp*erri + ki*mpii;
        
        Pchei = 1/COP*Pchi;
        
        Pfi = alpha_f*mi; % ???
        Pei = Pchei + Pfi + Forecast.Demand.E(i)*1000 + Pbat_charge_i - Pfcfi - Pbat_discharge_i;
        J1(i,1) = scaleCost(i,1)*Pei/1000 + scaleCost(i,3)*Pfcfi/eta_fce/1000;% - cas(i)*R_asi;

    end
    
    minimize ( sum(abs(J1)) )
    
    subject to
        
        for i = 1:nS
            
            hour = ceil(i*Plant.optimoptions.Resolution);

            if i == 1
                Ti = T0;
                Twi = Tw0;
                mpii = mpi0;
                Pchi = Pch0;
                Pfcfi = Pfcf0;
                Ebati = Ebat0;
            else
                Ti = x((i-2)*6+1);
                Twi = x((i-2)*6+2);
                mpii = x((i-2)*6+3);
                Pchi = x((i-2)*6+4);
                Pfcfi = x((i-2)*6+5);
                Ebati = x((i-2)*6+6);
            end

            T_refi = u((hour-1)*4+1);
            Pfcf_refi = u((hour-1)*4+2);
            Pbat_discharge_i = u((hour-1)*4+3);
            Pbat_charge_i = u((hour-1)*4+4);
            
            Tiplus = x((i-1)*6+1);
            Twiplus = x((i-1)*6+2);
%             mpiiplus = x((i-1)*6+3);
            Pchiplus = x((i-1)*6+4);
            Pfcfiplus = x((i-1)*6+5);
            Ebatiplus = x((i-1)*6+6);
            
            
            erri = T_refi - Ti;
            x((i-1)*6+3) == mpii + dt*erri;
            
%             R_asi = R_as(hour);

            mi = kp*erri + ki*mpii;
            Tmixi = alpha*Forecast.Temperature(i) + (1-alpha)*T_approx;
            Pcci = mi*Cp*(Tmixi-Ts);
            
            T_new = Ti + dt*1/Cr*(1/R*(Twi-Ti) + mi*Cp*(Ts-T_approx) + Forecast.Demand.H(i)*1000);
            Tw_new = Twi + dt*1/Cw*(1/R*(Ti-Twi) + 1/R*(Forecast.Temperature(i)-Twi));
%             mpi_new = mpii + dt*erri;
            Pch_new = Pchi + dt*1/tau_ch*(-Pchi + Pcci);
            Pfcf_new = Pfcfi + dt*1/tau_fc*(-Pfcfi + Pfcf_refi);
            Ebat_new = Ebati - dt*1/eta_d*Pbat_discharge_i + dt*eta_c*Pbat_charge_i;
            
            T_new == Tiplus;
            Tw_new == Twiplus;
%             mpi_new == mpiiplus;
            Pch_new == Pchiplus;
            Pfcf_new == Pfcfiplus;
            Ebat_new == Ebatiplus;
            
            mi >= 1;
            
%             0 <= R_asi <= 1e3;
%             R_as == 0;
            
%             T_low + kas*R_asi <= Tiplus <= T_high - kas*R_asi;
            T_low <= Tiplus <= T_high;
            0 <= Pfcfiplus <= Pfcf_max;
            0 <= Ebatiplus <= Ebat_max;

            Tref_low <= T_refi <= Tref_high;
            0 <= Pfcf_refi <= Pfcf_max;
            0 <= Pbat_discharge_i <= Pbat_max;
            0 <= Pbat_charge_i <= Pbat_max;
                        
        end
        
        Ebat_new == Ebat0;

cvx_end

GenDisp = zeros(nS+1,5);
GenDisp(1,:) = IC;
for i = 1:nS
    hour = ceil(i*Plant.optimoptions.Resolution);
    T(i,1) = x((i-1)*6+1); %room temperature in C
    Tw(i,1) = x((i-1)*6+2); %wall temperature in C
    mpi(i,1) = x((i-1)*6+3); %Integral term for PI controller for mass flow of HVAC in kg/s
    T_refi(i,1) = u((hour-1)*4+1);%temperature setpoint
    Pbat_discharge(i,1) = u((hour-1)*4+3);
    Pbat_charge(i,1) = u((hour-1)*4+4);
    erri = T_refi(i,1) - T(i,1);
    mi(i,1) = kp*erri + ki*mpi(i,1);%Mass flow in HVAC system
    
    GenDisp(i+1,3) = x((i-1)*6+5)/1000; %fuel cell electrical output
    GenDisp(i+1,4) = x((i-1)*6+4)/1000;%chiller output in kW of cooling
    GenDisp(i+1,5) = x((i-1)*6+6)/3600/1000;%battery state of charge in kWh
    GenDisp(i+1,1) = GenDisp(i,4)/COP + alpha_f*mi(i,1)/1000 + Forecast.Demand.E(i) + Pbat_charge(i,1)/1000 - GenDisp(i,3) - Pbat_discharge(i,1)/1000; %power to/from the electric grid in kW
%     R_as_vec(i+1,1) = R_as(hour); %ancillary service in kW 
end

end%Ends function NRELoptimization