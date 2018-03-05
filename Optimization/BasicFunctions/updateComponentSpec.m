function gen = updateComponentSpec(gen,param,value)
% UB is the new upper limit
% LB is the new lower limit
% dX_dt is the new ramp rate (kW/hr)
% Output is optional and can include new efficiency curves
UB = gen.Size;
if strcmp(gen.Type,'CHP Generator') || strcmp(gen.Type,'Electric Generator') || strcmp(gen.Type,'Hydrogen Generator')
    if isfield(gen.VariableStruct.Startup,'Electricity')
        LB = gen.VariableStruct.Startup.Electricity(end);
    elseif isfield(gen.VariableStruct.Startup,'DirectCurrent')
        LB = gen.VariableStruct.Startup.DirectCurrent(end);
    end
elseif strcmp(gen.Type,'Electrolyzer')
    LB = gen.VariableStruct.Startup.Hydrogen(end);
elseif strcmp(gen.Type,'Heater')
    LB = gen.VariableStruct.Startup.Heat(end);
elseif strcmp(gen.Type,'Chiller')
    LB = gen.VariableStruct.Startup.Cooling(end);
end
dX_dt = gen.VariableStruct.dX_dt;
p = eig(gen.VariableStruct.StateSpace.A);
w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
zeta = -real(p(1)+p(2))/(2*w0);
switch param
    case 'UB'
        UB = value;
    case 'LB'
        LB = value;
    case 'dX_dt'
        dX_dt = value;
    case 'Output'
        gen.Output = value;
    case 'w0'
        w0 = value;
    case 'zeta'
        zeta = value;
end
scale = UB/gen.Size;
LB = LB*scale;
dX_dt = dX_dt*scale;
T_peak = (UB-LB)/dX_dt*3600;

if strcmp(param,'w0')
    dX_dtNew = secondOrderResponse(gen,[]);
elseif zeta == 1%critical damping coefficient
        w0 = 1.5895*pi/(T_peak);
        gen.VariableStruct.StateSpace.A = [0 1; -(w0^2) -2*zeta*w0;];
        gen.VariableStruct.StateSpace.B = [0; (w0^2);];
        gen.VariableStruct.StateSpace.C = [1 0];
        gen.VariableStruct.StateSpace.D = 0;
        dX_dtNew = secondOrderResponse(gen,[]);
else
%% find natural frequency which achieves desired dX_dt
    error = 1;
    last_error = 1;
    while abs(error)>1e-5
        gen.VariableStruct.StateSpace.A = [0 1; -(w0^2) -2*zeta*w0;];
        gen.VariableStruct.StateSpace.B = [0; (w0^2);];
        gen.VariableStruct.StateSpace.C = [1 0];
        gen.VariableStruct.StateSpace.D = 0;
        dX_dtNew = secondOrderResponse(gen,[]);
        error = max(-.5,(dX_dt-dX_dtNew)/dX_dt);
        if sign(error) ~= sign(last_error)
            w0 = w0*(1+.5*error);
        else
            w0 = w0*(1+error);
        end
        last_error = error;
    end
end
gen.VariableStruct.dX_dt = dX_dtNew;
if isfield(gen.VariableStruct,'StartCost')
    gen.VariableStruct.StartCost = gen.VariableStruct.StartCost*scale;
end
F = fieldnames(gen.Output);
F = F(~strcmp(F,'Capacity'));
if any(ismember(F,'Electricity')) && any(gen.Output.Electricity>0)  && any(ismember(F,'Heat')) && any(gen.Output.Heat>0) 
    chp = true;
    F = {'Electricity'};
else
    chp = false;
end
gen.Size = UB;
gen.VariableStruct.Startup = [];
gen.VariableStruct.Shutdown = [];
for j = 1:1:length(F)
    gen.VariableStruct.Startup.Time = [0,1e3];
    gen.VariableStruct.Shutdown.Time = [0,1e3];
    if any(gen.Output.(F{j})>0)
        gen.VariableStruct.Startup.(F{j}) = [0,LB];
        gen.VariableStruct.Shutdown.(F{j}) = [LB,0];
        if ~isfield(gen.VariableStruct.Shutdown,'Input')
            input = LB/interp1(gen.Output.Capacity,gen.Output.(F{j}),LB/UB);
            gen.VariableStruct.Startup.Input = [0,input];
            gen.VariableStruct.Shutdown.Input = [input,0];
        end
        if chp
            heat = input*interp1(gen.Output.Capacity,gen.Output.Heat,LB/UB);
            gen.VariableStruct.Startup.Heat = [0,heat];
            gen.VariableStruct.Shutdown.Heat = [heat,0];
        end
    end
end