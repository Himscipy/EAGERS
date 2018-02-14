function dX_dt = secondOrderResponse(Gen,handles)
UB = Gen.Size;
if strcmp(Gen.Type,'CHP Generator') 
    Names = {'Electricity';};
%     Names = {'Electricity';'Heat';};
    LB = Gen.VariableStruct.Startup.Electricity(end);
elseif strcmp(Gen.Type,'Electric Generator')
    Names = {'Electricity';};
    LB = Gen.VariableStruct.Startup.Electricity(end);
elseif strcmp(Gen.Type,'Heater')
    Names = {'Heat'};
    LB = Gen.VariableStruct.Startup.Heat(end);
elseif strcmp(Gen.Type,'Chiller')
    Names = {'Cooling'};
    LB = Gen.VariableStruct.Startup.Cooling(end);
end

T_peak = (UB-LB)/Gen.VariableStruct.dX_dt*3600;
if T_peak>4*60
    Dt = 60;
elseif T_peak>3.6e4
    Dt = 3600;
else
    Dt = 1;
end
%convert continuous state-space to discrete time state-space
A = expm(Gen.VariableStruct.StateSpace.A.*Dt);
B = (Gen.VariableStruct.StateSpace.A)\(expm(Gen.VariableStruct.StateSpace.A.*Dt)-eye(2))*Gen.VariableStruct.StateSpace.B;
C = Gen.VariableStruct.StateSpace.C;
D = Gen.VariableStruct.StateSpace.D;
Time = 5*T_peak;
dX_dt = [];
while isempty(dX_dt)
    nS = round(Time/Dt)+1;
    u = UB*linspace(1,1,nS);
    [n,n2] = size(Gen.VariableStruct.StateSpace.C);
    X0 = zeros(n2,1);
    x0 = LB;
    for i = 1:1:n
        X0(find(Gen.VariableStruct.StateSpace.C(i,:),1,'first'))=x0(i);
    end
    [y,t] = ssSimulateDiscrete(A,B,C,D,X0,u,Dt);
    if any((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1)))) && ~any(abs((y(end-9:end,1)-y(end,1))/y(end,1))>1e-3)
        nR = min(nonzeros((1:1:nS+1)'.*((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))))));
        tRise = interp1(y(nR-1:nR,1)-y(1,1),t(nR-1:nR),.95*(u(1)-y(1,1)))/3600; %rise time in hours
        dX_dt = ((UB-LB).*(0.95)./tRise);
    else
        Time = Time*5;
    end
end
if ~isempty(handles)
    h1 = handles.ResponseRate;
    cla(h1)
    if Dt==1
        plot(h1,t,y);
        xlabel(h1,'Time (s)')
    elseif Dt==60
        plot(h1,t/60,y);
        xlabel(h1,'Time (min)')
    elseif Dt==3600
        plot(h1,t/3600,y);
        xlabel(h1,'Time (hr)')
    end
    ylabel(h1,'Power (kW)')
    % set(get(h1,'Ylabel'),'String','Power (kW)')
    legend(h1, Names)
end
end%Ends function secondOrderResponse

function [Y,T] = ssSimulateDiscrete(A,B,C,D,X0,u,Dt)
n = length(u);
T = (0:Dt:n*Dt)';
X = zeros(length(X0),n+1);
X(:,1) = X0;
Y = zeros(n+1,1);
Y(1) = C*X0;
for t = 1:1:n
    X(:,t+1) = A*X(:,t) + B*u(t);
    Y(t+1) = C*X(:,t+1) + D*u(t);
end
end%Ends function solverFixedStep
