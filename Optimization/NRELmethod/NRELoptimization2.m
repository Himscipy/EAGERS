function [GenDisp,OptSchedule] = NRELoptimization2(IC,T0,Forecast,scaleCost)
global Plant
% Specifically set up for a plant with 1 electric utility, 1 FC, 1 chiller 1 battery in that order
dt = Plant.optimoptions.Resolution*3600;% time step for discrete model, seconds
nS = Plant.optimoptions.Horizon/Plant.optimoptions.Resolution;% # of time steps in optimization

c_as = 0.1*ones(nS,1); % ancillary service price
alpha = 0.4; % outdoor air ratio
Cp_air = 1e3; % specific heat of air
Ts = 12.778; % supply air temperature
T_approx = 22.22; % approximation for room temperature for complexity
COP = 3.5; 
tau_fc = 120*60; % time constant of the FC, second
alpha_f = 500; % fan power coefficient
kas = 1e-4; % coefficent to calculate temperature buffer from ancillary service

T_high = 23.33;
T_low = 21.11;

m_min = 1;
m_max = 10;

Pfc_f0 = IC(3)*1000/Plant.Generator(3).Output.Electricity(end);%initial fuel flow of the fuel cell in Watts
Eb0 = IC(5)*3.6e6;%stored energy in battery in Joules
%%

cvx_begin
% cvx_solver_settings( 'dumpfile', 'test' ) 
cvx_solver_settings( 'NumericFocus', '3' ) 
% cvx_precision best

    % Decision variables
    variables T_ref(nS) Pb_in(nS) Pb_out(nS) m(nS) Pfc_f_ref(nS) Pfc_h(nS) Pboil(nS) Ras(nS)
    % State variables
    variables T(nS) Eb(nS) Pfc_f(nS)
    
    for i = 1:nS
        Tmix_i = alpha*Forecast.Weather.Tdb(i) + (1-alpha)*T_approx;
        Pgrid = (m(i)*Cp_air*(Tmix_i-Ts)/COP + alpha_f*m(i) + Pb_in(i) + Forecast.Demand.E(i)*1000 - Pfc_f(i)*Plant.Generator(3).Output.Electricity(end) - Pb_out(i));%Grid power in Watts
        J(i,1) = scaleCost(i,1)*Pgrid + scaleCost(i,2)*(Pfc_f(i) + Pboil(i)/Plant.Generator(6).Output.Heat(end)) - c_as(i)*Ras(i);%Cost for this generation & fuel
    end

    minimize ( sum(abs(J)) )
    
    subject to
        
        for i = 1:nS
            if i == 1
                T_i = T0;
                Pfc_f_i = Pfc_f0;
                Eb_i = Eb0;
            else
                T_i = T(i-1);
                Pfc_f_i = Pfc_f(i-1);
                Eb_i = Eb(i-1);
            end
            T_iplus = T(i);
            Pfc_f_iplus = Pfc_f(i);
            Eb_iplus = Eb(i);

            T_ref_i = T_ref(i);
            Pb_in_i = Pb_in(i);
            Pb_out_i = Pb_out(i);
            m_i = m(i);
            Pfc_f_ref_i = Pfc_f_ref(i);
            Pfc_h_i = Pfc_h(i);
            Pboil_i = Pboil(i);
            Ras_i = Ras(i);
            
            Qcool_i = m_i*Cp_air*(Ts-T_approx); % HVAC cooling into the room
            Qheat_i = Pfc_h_i + Pboil_i; % HVAC heating into the room
            
            T_new = 0.99231*T_i + 1.0556e-11*Forecast.Weather.Tdb(i) + 4.9938e-6*(-Qcool_i+Qheat_i+Forecast.Demand.IntGain(i)) + 2.1844e-4*Forecast.Weather.irradDireNorm(i);
            Pfc_f_new = Pfc_f_i + dt*1/tau_fc*(-Pfc_f_i + Pfc_f_ref_i);
            Eb_new = Eb_i + dt*Plant.Generator(5).VariableStruct.ChargeEff*Pb_in_i - dt*1/Plant.Generator(5).VariableStruct.DischargeEff*Pb_out_i;
            
            T_iplus == T_new;
            Pfc_f_iplus == Pfc_f_new;
            Eb_iplus == Eb_new;
            
            T_iplus == T_ref_i;
            
            % Limits
            Pboil_i >= 0;
            0 <= Pfc_h_i <= Pfc_f_i*Plant.Generator(3).Output.Heat(end);
            0 <= Ras_i <= 1e3;
            m_min <= m_i <= m_max;
            T_low + kas*Ras_i <= T_iplus <= T_high - kas*Ras_i;
            0 <= Pfc_f_iplus <= Plant.Generator(3).Size/Plant.Generator(3).Output.Electricity(end)*1000;
            0 <= Eb_iplus <= Plant.Generator(5).Size*3.6e6;

            -18 <= T_ref_i <= 38;
            0 <= Pfc_f_ref_i <= Plant.Generator(3).Size/Plant.Generator(3).Output.Electricity(end)*1000;
            0 <= Pb_in_i <= Plant.Generator(5).VariableStruct.PeakCharge*3.6e6;
            0 <= Pb_out_i <= Plant.Generator(5).VariableStruct.PeakDisch*3.6e6;
        
        end
        
        Eb_new == Eb0;

cvx_end

GenDisp = zeros(nS+1,6);
GenDisp(1,:) = IC;
for i = 1:nS
    Tmix_i = alpha*Forecast.Weather.Tdb(i) + (1-alpha)*T_approx;
    GenDisp(i+1,3) = Pfc_f(i,1)*Plant.Generator(3).Output.Electricity(end)/1000; %fuel cell electrical output
    GenDisp(i+1,4) = (m(i)*Cp_air*(Tmix_i-Ts))/1000;%chiller output in kW of cooling
    GenDisp(i+1,5) = GenDisp(i,5) + (dt*Plant.Generator(5).VariableStruct.ChargeEff*Pb_in(i,1) - dt*1/Plant.Generator(5).VariableStruct.DischargeEff*Pb_out(i))/3.6e6;%battery state of charge in kWh
    GenDisp(i+1,6) = Pboil(i)/1000;
    GenDisp(i+1,1) = GenDisp(i,4)/COP + alpha_f*m(i,1)/1000 + Forecast.Demand.E(i) + Pb_in(i,1)/1000 - GenDisp(i,3) - Pb_out(i,1)/1000; %power to/from the electric grid in kW
end
OptSchedule.T_ref = T_ref;
OptSchedule.Pfc_h = Pfc_h;
