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
SS = ss(Gen.VariableStruct.StateSpace.A,Gen.VariableStruct.StateSpace.B,Gen.VariableStruct.StateSpace.C,Gen.VariableStruct.StateSpace.D);
SS = c2d(SS,Dt);
Time = 5*T_peak;
dX_dt = [];
while isempty(dX_dt)
    nS = round(Time/Dt)+1;
    t = linspace(0, Dt*(nS-1),nS);
    u = UB*linspace(1,1,nS);
    [n,n2] = size(Gen.VariableStruct.StateSpace.C);
    X0 = zeros(n2,1);
    x0 = LB;
    for i = 1:1:n
        X0(find(Gen.VariableStruct.StateSpace.C(i,:),1,'first'))=x0(i);
    end
    [y,t] = lsim(SS,u,t,X0);
    if any((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1)))) && ~any(abs((y(end-9:end,1)-y(end,1))/y(end,1))>1e-3)
        nR = min(nonzeros((1:1:nS)'.*((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))))));
        tRise = interp1(y(nR-1:nR,1)-y(1,1),t(nR-1:nR),.95*(u(1)-y(1,1)))/3600; %rise time in hours
%         tRise = t(find((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))),1,'first'))/3600; %rise time in hours
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
