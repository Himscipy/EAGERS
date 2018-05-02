function GenDisp = NRELoptimization(IC,Forecast,scaleCost)
global Plant

dt = Plant.optimoptions.Resolution*3600;% time step for discrete model, seconds
nS = Plant.optimoptions.Horizon/Plant.optimoptions.Resolution;% # of time steps in optimization
N_schedule = Plant.optimoptions.Horizon; % assumes hourly scheduling
Forecast.Demand.H = Forecast.Demand.H;


InitialStates = [22;0;0;0.5*1e3*3600;0];
tshift = 24*10;

% Initial conditions for building and equipment status
T0 = InitialStates(1);
Pfc_f0 = IC(2);
Pfc_h0 = InitialStates(3);
Eb0 = IC(5);
u0 = InitialStates(5);

%dt = 3*60; % Time step for discrete thermal dynamics model, seconds
dt_h = 60*60; % Time step for scheduling period
%N_sch = 24; % Number of scheduling period
N_state = 3600/dt*N_schedule; % Number of thermal dynamics model steps
N_perhour = dt_h/dt; % Number of thermal dynamics model steps per hour


Toa = Forecast.Temperature(1+tshift:N_sch+tshift); % outside air temperature
Pl_h_6min = Forecast.Demand.H*1e3;
tvec_6min = (6:6:8760*60)';
tvec_1hour = (60:60:8760*60)';
Pl_h_1hour = interp1(tvec_6min,Pl_h_6min,tvec_1hour);

Pl_h = Pl_h_1hour(1+tshift:N_sch+tshift);
Pl_e = Forecast.Demand.E(1+tshift:N_sch+tshift);

% Prices
c_e_perday = 0.1*ones(N_sch,1)*1e-3; % electricity price
c_e_perday(N_sch*2/3:N_sch*5/6) = 0.2*1e-3;
c_e_120day = repmat(c_e_perday,365,1);
c_e = c_e_120day(1+tshift:N_sch+tshift);
c_ng_perday = 0.03*ones(N_sch,1)*1e-3; % natural gas price
c_ng_perday(1:N_sch/3) = 0.03*1e-3;
c_ng_perday(N_sch*5/6:end) = 0.03*1e-3;
c_ng_120day = repmat(c_ng_perday,365,1);
c_ng = c_ng_120day(1+tshift:N_sch+tshift);
c_as_perday = 0.03*ones(N_sch,1)*1e-3; % ancillary service price
c_as_120day = repmat(c_as_perday,365,1);
c_as = c_as_120day(1+tshift:N_sch+tshift);

% alpha = 0.4; % outdoor air ratio
m_oa = 0.3; % outdoor air flow rate, kg/s


% % Building zone parameters
Zoneparam.Cp = 1e3; % specific heat of air
Zoneparam.Cr = 5e6;
Zoneparam.Cw = 3.75e8;
Zoneparam.R1 = 0.031e-3;
Zoneparam.R2 = 0.5898e-3;
Zoneparam.Ts = f2c(55); % supply air temperature
Zoneparam.T_approx = 22; % approximation for room temperature for complexity

COP = Plant.Generator(4).Output.Cooling(end); % chiller COP
eta_fc_e = Plant.Generator(3).Output.Electricity(end); % electrical efficiency of the FC
Pfc_f_max = Plant.Generator(3).Size*1000; % FC fuel consumption capacity
tau_ch = 200; % time constant of the chiller, second
tau_fc = 120*60; % time constant of the FC, second


DischCurrent = Plant.Generator(5).VariableStruct.PeakDisch.*Plant.Generator(5).Size/Plant.Generator(5).VariableStruct.Voltage*1000;
DischResistScaled = (100/DischCurrent)*Plant.Generator(5).VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
DischVoltLoss = DischCurrent*DischResistScaled; %keep in mind when calculating loss as function of discharge current
ChargeCurrent = Plant.Generator(5).VariableStruct.PeakCharge*Plant.Generator(5).Size/Plant.Generator(5).VariableStruct.Voltage*1000;
ChargeResistScaled = (100/ChargeCurrent)*Plant.Generator(5).VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
ChargeVoltLoss = ChargeCurrent*ChargeResistScaled;
eta_c = Plant.Generator(5).VariableStruct.Voltage/(Plant.Generator(5).VariableStruct.Voltage+ChargeVoltLoss); %charging efficiency
eta_d = (Plant.Generator(5).VariableStruct.Voltage-DischVoltLoss)/Plant.Generator(5).VariableStruct.Voltage; %discharging efficiency
Eb_max = Plant.Generator(5).Size*(Plant.Generator(5).VariableStruct.MaxDOD/100)*1000*3600; % usable battery capacity (J)
Pb_max = (DischCurrent*Plant.Generator(5).VariableStruct.Voltage*eta_d); % battery ramping capacity (J/s)
eta_fc_h = 0.35; % heat efficiecny of the FC
eta_boil = 0.9; % boiler efficiency
a_f = 2.225e+03; % fan power coefficient
b_f = -2.977e+03;

kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service
c_up = 1; % FC start-up cost
c_down = 1; % FC shut-down cost



T_high = 23; % temperature bounds
T_low = 21;
Tref_high = 32; % temperature set-points bounds
Tref_low = 12;
J = zeros(nS,1);

Pfc_f_min = 0.2*Pfc_f_max; % FC fuel consumption minimum operating point

m_min = 2.23; % supply air flow rate limits
m_max = 10;
FC_ramp_up = 0.25*Pfc_f_max; % FC ramp rate limits
FC_ramp_down = 0.25*Pfc_f_max;
min_up = 6; % FC minimum up time
min_down = 6; % FC minimum down time

[num_vec,den_vec]=ZoneThermalModelEDC(Zoneparam,dt); % Model parameter generator

%% EDC optimization
cvx_begin

% Desicion variables:
% zone temperatue set-point, battery charge, battery discharge, supply
% air flow rate, FC fuel consumption set-point, FC heat output
% set-point, total reheat, AS capacity
%variables T_ref(N_sch) Pb_in(N_sch) Pb_out(N_sch) m(N_sch) Pfc_f_ref(N_sch) Pfc_h_ref(N_sch) Prh(N_sch) Ras(N_sch)
variables u(8*N_schedule)
% State variables:
% zone temperature, battery SOC
%variables T(N_state) Eb(N_state)
variables x(2*N_state)
% Interger variables:
% FC on/off, FC start up, FC shut down
variable w(3*(N_schedule+max(min_up,min_down))) binary

% calculate cost
for i = 1:nS
    hour = ceil(i*Plant.optimoptions.Resolution);
    
 
    Pb_in_i = u((hour-1)*8+2);
    Pb_out_i = u((hour-1)*8+3);
    m_i = u((hour-1)*8+4);
    Pfc_f_ref_i = u((hour-1)*8+5);
    Pfc_h_ref_i = u((hour-1)*8+6);
    Prh_i = u((hour-1)*8+7);
    Ras_i = u((hour-1)*8+8);
    
    u_up_i = w((hour-1)*3+2);
    u_down_i = w((hour-1)*3+3);
    
% Total cost incurred from various sources (electricity, natural gas, fuel cell, AS)
    % Electricity cost
    Pcc_i = (m_i-m_oa)*Zoneparam.Cp*(Zoneparam.T_approx-Zoneparam.Ts) + m_oa*Cp*(Toa(i)-Zoneparam.Ts); % cooling coil heat exchange
    Pch_i = 1/COP*Pcc_i; % chiller power
    
    Pf_i = a_f*m_i + b_f; % fan power
    
    Pfc_e_i = Pfc_f_ref_i*eta_fc_e; % FC elec output
    
    P_e_i = Pch_i + Pf_i + Pb_in_i + Forecast.Demand.E(i)- Pfc_e_i - Pb_out_i; % grid power
    
    J_e_i = scaleCost(i,1)*P_e_i;
    % Gas cost
    P_ng_i = Pfc_f_ref_i + (Prh_i-Pfc_h_ref_i)/eta_boil;
    J_gas_i = scaleCost(i,2)*P_ng_i;    
    % AS payment
    J_as_i = -scaleCost(i,3)*Ras_i;    
    % FC on/off
    J_fc_com_i = c_up*u_up_i + c_down*u_down_i;   
    % Total cost
    J(i,1) = J_e_i + J_gas_i + J_as_i + J_fc_com_i;
    
end


minimize ( sum(abs(J)) )

subject to

for i_sch = 1:nS  % (nS=24)
    
    hour = ceil(i*Plant.optimoptions.Resolution);
    
    T_ref_i = u((hour-1)*8+1);
    Pb_in_i = u((hour-1)*8+2);
    Pb_out_i = u((hour-1)*8+3);
    m_i = u((hour-1)*8+4);
    Pfc_f_ref_i = u((hour-1)*8+5);
    Pfc_h_ref_i = u((hour-1)*8+6);
    Prh_i = u((hour-1)*8+7);
    Ras_i = u((hour-1)*8+8);
    u_i =  w((hour-1)*3+1);
    
    if i_sch == 1
        Pfc_f_ref_iminus = Pfc_f0;
        Pfc_h_ref_iminus = Pfc_h0;
        u_iminus = u0;
    else
        Pfc_f_ref_iminus = u((hour-1)*8+5-8);
        Pfc_h_ref_iminus = u((hour-1)*8+6-8);
       u_iminus = w((hour-1)*3+1-3);
    end
    
    % Limits
    0 <= Pfc_h_ref_i <= Pfc_f_ref_i*eta_fc_h;
    0 <= Prh_i - Pfc_h_ref_i;
    0 <= Ras_i <= 1e3;
    m_min <= m_i <= m_max;
    Pfc_f_min*u_i <= Pfc_f_ref_i <= Pfc_f_max*u_i;
    Pfc_f_min*1 <= Pfc_f_ref_i <= Pfc_f_max*1;
    -FC_ramp_down <= Pfc_f_ref_i - Pfc_f_ref_iminus <= FC_ramp_up;
    0 <= Pb_in_i <= Pb_max;
    0 <= Pb_out_i <= Pb_max;
    Tref_low <= T_ref_i <= Tref_high;
    
    % FC commitment
    for i_up = 1:min_up % FC minimum up time
        -u_iminus + u_i - w((hour-1)*3+1+3*i_up) <= 0;
    end
    for i_down = 1:min_down % FC minimum down time
        -u_iminus + u_i - w((hour-1)*3+1+3*i_down) <= 0;
    end
    
    -u_iminus + u_i - w((hour-1)*3+2) <= 0; % FC start up
    u_iminus - u_i - w((hour-1)*3+3) <= 0; % FC shut down
%     
    % States in faster time step
    for i_state = 1:N_perhour
        T_i = x((i_state-1)*2+1+N_perhour*2*(i_sch-1));
        Eb_i = x((i_state-1)*2+2+N_perhour*2*(i_sch-1));
        if i_sch == 1 && i_state == 1
            T_iminus = T0;
            T_iminus2 = T0;
            %                     Tw_iminus = Tw0;
            Eb_iminus = Eb0;
        elseif i_sch == 1 && i_state == 2
            T_iminus = x((i_state-2)*2+1+N_perhour*2*(i_sch-1));
            T_iminus2 = T0;
            %                     Tw_iminus = Tw0;
            Eb_iminus = x((i_state-1)*2+N_perhour*2*(i_sch-1));
        else
       
            T_iminus = x((i_state-1)*2-1+N_perhour*2*(i_sch-1));
            T_iminus2 = x((i_state-1)*2-3+N_perhour*2*(i_sch-1));
            Eb_iminus = x((i_state-1)*2+N_perhour*2*(i_sch-1));
        end
        Pl_h_iminus = Pl_h(i_sch);
        Pl_h_iminus2 = Pl_h(i_sch);
        Toa_iminus = Toa(i_sch);
        Toa_iminus2 = Toa(i_sch);
        
        
        T_new = -den_vec(3,2)*T_iminus - den_vec(3,3)*T_iminus2 + num_vec(1,2)*m_i + num_vec(2,2)*(Pl_h_iminus+Prh_i) + num_vec(3,2)*Toa_iminus + num_vec(1,3)*m_i + num_vec(2,3)*(Pl_h_iminus2+Prh_i) + num_vec(3,3)*Toa_iminus2;
        Eb_new = Eb_iminus + dt*eta_c*Pb_in_i - dt*1/eta_d*Pb_out_i;
        
        T_new == T_i;
                         
        Eb_new == Eb_i
        
        T_low + kas*Ras_i <= T_i <= T_high - kas*Ras_i;
        0 <= Eb_i <= Eb_max;
        
    end
    
    T_new == T_ref_i; % temperature reach set-point
    
    
end

Eb_new == Eb0; % final battery SOC

cvx_end

GenDisp = zeros(nS,4); 
for i = 1:nS
    hour = ceil(i*Plant.optimoptions.Resolution);
    Pb_in_i(i,1) = u((hour-1)*8+2);
    Pb_out_i(i,1) = u((hour-1)*8+3);
    m_i(i,1) = u((hour-1)*8+4);
    Pfc_f_ref_i(i,1) = u((hour-1)*8+5);
    Pfc_h_ref_i(i,1) = u((hour-1)*8+6);
    Prh_i(i,1) = u((hour-1)*8+7);
    Ras_i(i,1) = u((hour-1)*8+8);
    
    Pcc_i(i,1) = (m_i(i,1)-m_oa)*Zoneparam.Cp*(Zoneparam.T_approx-Zoneparam.Ts) + m_oa*Cp*(Toa(i)-Zoneparam.Ts); % cooling coil heat exchange
    Pch_i(i,1) = 1/COP*Pcc_i(i,1);
    
    Pf_i(i,1) = a_f*m_i(i,1) + b_f;
    
    Pfc_e_i(i,1) = Pfc_f_ref_i(i,1)*eta_fc_e;
    
    u_up_i(i,1) = w((hour-1)*3+2);
    u_down_i(i,1) = w((hour-1)*3+3);
    GenDisp(i+1,2)= Pfc_e_i(i,1); % Fuel cell electrical output
    GenDisp(i+1,3)= Pcc_i(i,1);   % Chiller output
    GenDisp(i+1,4)= x(2+N_perhour*2*(i-1));   %battery state of charge
    GenDisp(i+1,1) = Pch_i(i,1) + Pf_i(i,1) + Pb_in_i(i,1) + Forecast.Demand.E(i)- Pfc_e_i(i,1) - Pb_out_i(i,1)+...
    Pfc_f_ref_i(i,1) + (Prh_i(i,1)-Pfc_h_ref_i(i,1))/eta_boil-Ras_i(i,1)+u_up_i(i,1) + u_down_i(i,1);%power to/from the electric grid in kW  
end







