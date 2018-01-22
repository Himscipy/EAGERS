function component = updateComponentSpec(component,dX_dt_Set,UB,LB,Output)
%load the component you wish to modify into component
% i is the index in Plant.Generator that you will put the updated component
% UB is the new upper limit
% LB is the new lower limit
% dX_dt is the new ramp rate (kW/hr)
% Output is optional and can include new efficiency curves

% global Plant
T_peak = (UB-LB)/dX_dt_Set*3600;
w0 = 1.5895*pi/(T_peak);
zeta =1; %damping coefficient
A = [0 1; -(w0^2) -2*zeta*w0;];
B = [0; (w0^2);];
C = [1 0];
D = 0;

component.VariableStruct.StateSpace.A = A;
component.VariableStruct.StateSpace.B = B;
component.VariableStruct.StateSpace.C = C;
component.VariableStruct.StateSpace.D = D;
%remove field Dt
if isfield(component.VariableStruct.StateSpace, 'Dt')
    component.VariableStruct.StateSpace = rmfield(component.VariableStruct.StateSpace,'Dt');
end

if isempty(Output)
    F = fieldnames(component.Output);
    F = F(~strcmp(F,'Capacity'));
    Output.Capacity = component.Output.Capacity;
    for j = 1:1:length(F)
        if any(component.Output.(F{j})>0)
            Output.(F{j}) = component.Output.(F{j});
        end
    end
end

F = fieldnames(Output);
F = F(~strcmp(F,'Capacity'));
if any(ismember(F,'Electricity')) && any(Output.Electricity>0)  && any(ismember(F,'Heat')) && any(Output.Heat>0) 
    chp = true;
    F = {'Electricity'};
else chp = false;
end
component.Size = UB;
component.VariableStruct.Startup = [];
component.VariableStruct.Shutdown = [];
for j = 1:1:length(F)
    component.VariableStruct.Startup.Time = [0,1e3];
    component.VariableStruct.Shutdown.Time = [0,1e3];
    if any(Output.(F{j})>0)
        component.VariableStruct.Startup.(F{j}) = [0,LB];
        component.VariableStruct.Shutdown.(F{j}) = [LB,0];
        if ~isfield(component.VariableStruct.Shutdown,'Input')
            input = LB/interp1(Output.Capacity,Output.(F{j}),LB/UB);
            component.VariableStruct.Startup.Input = [0,input];
            component.VariableStruct.Shutdown.Input = [input,0];
        end
        if chp
            heat = input*interp1(Output.Capacity,Output.Heat,LB/UB);
            component.VariableStruct.Startup.Heat = [0,heat];
            component.VariableStruct.Shutdown.Heat = [heat,0];
        end
    end
end
component.Output = Output;
% Plant.Generator(i) = component;

% %% Plot response
% Dt = 1;
% SS = ss(A,B,C,D);
% SS = c2d(SS,Dt);
% Time = 10*T_peak;
% nS = round(Time/Dt)+1;
% t = linspace(0, Dt*(nS-1),nS);
% u = UB*linspace(1,1,nS);
% [n,n2] = size(C);
% X0 = zeros(n2,1);
% x0 = LB;
% for i = 1:1:n
%     X0(find(C(i,:),1,'first'))=x0(i);
% end
% [y,t] = lsim(SS,u,t,X0);
% plot(t/60,y);
% xlabel('Time (min)')
% 
% tRise = t(find((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))),1,'first'))/3600; %rise time in hours
% dX_dt = ((UB-LB).*(0.95)./tRise)