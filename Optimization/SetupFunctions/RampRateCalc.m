function [dX_dt, SS_1] = RampRateCalc(SS,LB,UB,H_0)
SS = SS(1,:);
SS_1 = SS;
if isfield(SS,'Dt') && ~isempty(SS.Dt)
    Dt = SS.Dt;
    SS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
else
    SS_1.Dt = [];
    Dt =1;
    SS.Dt = 1;
    SS = ss(SS(1).A,SS(1).B,SS(1).C,SS(1).D);
    SS = c2d(SS,Dt);
end
if Dt~=1 %convert to 1 second sampling time
    SS = d2d(SS,1);
    r = length(SS);
    SS_1.A = SS(1).A;
    SS_1.B = SS(1).B;
    SS_1.C = SS(1).C;
    SS_1.D = SS(1).D;
    for k = 2:1:r
        SS_1.C(end+1,:) = SS(k).C;
        SS_1.D(end+1,:) = SS(k).D;
    end
end
[z,z2] = size(SS.C);
x0 = LB;
if ~isempty(H_0)
    x0(2) = H_0; %CHP heat produced per unit electricity
elseif z>1 %if there is no heat ratio, but there is multiple state outputs, then it is capable of CHP, but is not connected for CHP
    x0(2) = 0;
end
nS = round(12*3600/Dt)+1; % assume ramping is less than 4 hours (i forget why I made this limit)
t = linspace(0, Dt*(nS-1),nS);
u = UB*linspace(1,1,nS)';
X0 = zeros(z2,1);
for k = 1:1:z
    X0(find(SS.C(k,:),1,'first'))=x0(k);%finds the first non-zero element in SS.C and makes X0 at that index = x0
end
SS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
[y,t] = lsim(SS,u,t,X0);
tRise = t(find((y(:,1)-y(1,1))>(.95*(u(1)-y(1,1))),1,'first'))/3600; %rise time in hours
if isempty(tRise)
    tRise = 12;
end
dX_dt = ((UB-LB).*(0.95)./tRise);
