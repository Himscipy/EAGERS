function [GenDisp,Feasible] = DispatchQP(QP,Locked)
[m,n] = size(QP.organize);
nG = length(QP.constCost);
nS = m-1;

QP = disableGenerators(QP,Locked,[]);%Disable generators here
[GenSetting, cost,Feasible] = callQPsolver(QP);
GenDisp = zeros(nS+1,n);
if Feasible ~=1
    Feasible = false;
%     disp('Infeasible in DispatchQP');
else
    Feasible = true;
    for i = 1:1:nG
        if any(QP.Renewable(:,i)~=0)
            GenDisp(1,i) = RenewableOutput(i,[],'Actual');
            GenDisp(2:end,i) = QP.Renewable(:,i);
        else
            for t = 1:1:nS+1
                if ~isempty(QP.organize{t,i})
                    GenDisp(t,i) = sum(GenSetting(QP.organize{t,i}));%%put solution back into recognizable format
                end
            end
        end
    end
    for i = nG+1:1:n %Transmission lines
        if ~isempty(QP.organize{1,i})
            GenDisp(1,i) = GenSetting(QP.organize{1,i});%%put solution back into recognizable format
        end
        for t = 2:1:nS+1
            GenDisp(t,i) = GenSetting(QP.organize{t,i});%%put solution back into recognizable format
        end
    end
    GenDisp(abs(GenDisp)<1e-3) = 0;
end
% % Check charging state
% global Plant
% for i = 1:1:nG
%     n = 0;
%     if ~isempty(QP.Organize.Inequalities(i))
%         n = n+1;
%         states = eval(QP.Organize.States(i));
%         chargeState(:,n) = zeros(nS,1);
%         for t = 1:1:nS
%             chargeState(t,1) = GenSetting(QP.organize{t+1,i}+1);%%charge state (if it exists) is adjacent to SOC
%         end
%         eff = Plant.Generator(i).OpMatA.Stor.DischEff*Plant.Generator(i).OpMatA.Stor.ChargeEff;
%         StorPower(:,n) = GenDisp(1:nS,i)-GenDisp(2:nS+1,i);
%         Charge = -StorPower(:,n).*(StorPower(:,n)<0);
%         correctChargePower = Charge*(1-eff);
%         error(:,n) = chargeState(:,n)- correctChargePower;
%     end
% end