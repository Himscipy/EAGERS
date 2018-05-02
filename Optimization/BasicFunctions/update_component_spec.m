function gen = update_component_spec(gen,param,value)
% upper_bound is the new upper limit
% lower_bound is the new lower limit
% dx_dt is the new ramp rate (kW/hr)
upper_bound = gen.Size;
if strcmp(gen.Type,'CHP Generator') || strcmp(gen.Type,'Electric Generator') || strcmp(gen.Type,'Hydrogen Generator')
    if isfield(gen.VariableStruct.Startup,'Electricity')
        lower_bound = gen.VariableStruct.Startup.Electricity(end);
    elseif isfield(gen.VariableStruct.Startup,'DirectCurrent')
        lower_bound = gen.VariableStruct.Startup.DirectCurrent(end);
    end
elseif strcmp(gen.Type,'Electrolyzer')
    lower_bound = gen.VariableStruct.Startup.Hydrogen(end);
elseif strcmp(gen.Type,'Heater')
    lower_bound = gen.VariableStruct.Startup.Heat(end);
elseif strcmp(gen.Type,'Chiller')
    lower_bound = gen.VariableStruct.Startup.Cooling(end);
elseif strcmp(gen.Type,'Cooling Tower')
    lower_bound = gen.VariableStruct.Startup.heat_reject(end);
end
dx_dt = gen.VariableStruct.dX_dt;
p = eig(gen.VariableStruct.StateSpace.A);
w0 = sqrt(real(p(1))^2 + imag(p(1))^2);
zeta = -real(p(1)+p(2))/(2*w0);
switch param
    case 'UB'
        upper_bound = value;
    case 'LB'
        lower_bound = value;
    case 'dX_dt'
        dx_dt = value;
    case 'Output'
        gen.Output = value;
    case 'w0'
        w0 = value;
    case 'zeta'
        zeta = value;
end
scale = upper_bound/gen.Size;
lower_bound = lower_bound*scale;
dx_dt = dx_dt*scale;
t_peak = (upper_bound-lower_bound)/dx_dt*3600;

F = fieldnames(gen.Output);
F = F(~strcmp(F,'Capacity'));
if any(ismember(F,'Electricity')) && any(gen.Output.Electricity>0)  && any(ismember(F,'Heat')) && any(gen.Output.Heat>0) 
    chp = true;
    F = {'Electricity'};
else
    chp = false;
end
gen.Size = upper_bound;
gen.VariableStruct.Startup = [];
gen.VariableStruct.Shutdown = [];
for j = 1:1:length(F)
    gen.VariableStruct.Startup.Time = [0,1e3];
    gen.VariableStruct.Shutdown.Time = [0,1e3];
    if any(gen.Output.(F{j})>0)
        gen.VariableStruct.Startup.(F{j}) = [0,lower_bound];
        gen.VariableStruct.Shutdown.(F{j}) = [lower_bound,0];
        if ~isfield(gen.VariableStruct.Shutdown,'Input')
            input = lower_bound/interp1(gen.Output.Capacity,gen.Output.(F{j}),lower_bound/upper_bound);
            gen.VariableStruct.Startup.Input = [0,input];
            gen.VariableStruct.Shutdown.Input = [input,0];
        end
        if chp
            heat = input*interp1(gen.Output.Capacity,gen.Output.Heat,lower_bound/upper_bound);
            gen.VariableStruct.Startup.Heat = [0,heat];
            gen.VariableStruct.Shutdown.Heat = [heat,0];
        end
    end
end

if strcmp(param,'w0')
    dx_dt_new = ss_response(gen,[]);
elseif zeta == 1%critical damping coefficient
        w0 = 1.5895*pi/(t_peak);
        gen.VariableStruct.StateSpace.A = [0 1; -(w0^2) -2*zeta*w0;];
        gen.VariableStruct.StateSpace.B = [0; (w0^2);];
        gen.VariableStruct.StateSpace.C = [1 0];
        gen.VariableStruct.StateSpace.D = 0;
        dx_dt_new = ss_response(gen,[]);
else
%% find natural frequency which achieves desired dX_dt
    error = 1;
    last_error = 1;
    while abs(error)>1e-5
        gen.VariableStruct.StateSpace.A = [0 1; -(w0^2) -2*zeta*w0;];
        gen.VariableStruct.StateSpace.B = [0; (w0^2);];
        gen.VariableStruct.StateSpace.C = [1 0];
        gen.VariableStruct.StateSpace.D = 0;
        dx_dt_new = ss_response(gen,[]);
        error = max(-.5,(dx_dt-dx_dt_new)/dx_dt);
        if sign(error) ~= sign(last_error)
            w0 = w0*(1+.5*error);
        else
            w0 = w0*(1+error);
        end
        last_error = error;
    end
end
gen.VariableStruct.dX_dt = dx_dt_new;
if isfield(gen.VariableStruct,'StartCost')
    gen.VariableStruct.StartCost = gen.VariableStruct.StartCost*scale;
end