function [Dispatch,LineLoss,excessHeat,excessCool,hydroSOC,Temperature,Heating,Cooling] = sortSingleSolution(x,QP)
[~,n] = size(QP.organize);
nG = length(QP.constCost);
nB = length(QP.Organize.Building.r);
nL = n-nG-nB;
nH = nnz(QP.Organize.Hydro);
hydroSOC = zeros(1,nH);
excessHeat = zeros(1,nnz(QP.Organize.HeatVented));
excessCool = zeros(1,nnz(QP.Organize.CoolVented));
LineFlows = zeros(1,nL);
LineLoss = zeros(1,nL);
Temperature = zeros(1,nB);
Heating = zeros(1,nB);
Cooling = zeros(1,nB);
GenDisp = zeros(1,nG);
for i = 1:1:nG
    if isfield(QP,'Renewable') && any(QP.Renewable(:,i)~=0)
        GenDisp(2:end,i) = QP.Renewable(1,i);
    else
        Out_vs_State = QP.Organize.Out_vs_State{1,i};%linear maping between state and output
        if ~isempty(QP.organize{1,i})
            states = QP.organize{1,i};
            for j = 1:1:length(states)
                GenDisp(1,i) = GenDisp(1,i) + Out_vs_State(j)*x(states(j)); %record this combination of outputs (SOC for storage)
            end
        end
    end
    
end
for i = 1:1:length(QP.Organize.Hydro)
    hydroSOC(1,i) = x(QP.organize{1,QP.Organize.Hydro(i)}+1,1);
end
GenDisp(abs(GenDisp)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors

for i = 1:1:nL
    LineFlows(1,i) = sum(x(QP.organize{1,i+nG}));%power transfer
    if length(QP.Organize.States{i+nG})>1
        LineLoss(1,i) = sum(x(QP.organize{1,i+nG}+1)); %down (positive) lines
        LineLoss(1,i) = LineLoss(1,i) + sum(x(QP.organize{1,i+nG}+2)); %up (negative) lines
    end
end

for i = 1:1:nB
    Tset = QP.organize{1,i+nG+nL};
    Temperature(1,i) = x(Tset,1);
    Heating(1,i) = x(Tset+1,1) - QP.Organize.Building.H_Offset(1,i);
    Cooling(1,i) = x(Tset+2,1) - QP.Organize.Building.C_Offset(1,i);
end
Dispatch = [GenDisp,LineFlows,Temperature];
%pull out any dumped heat
for i = 1:1:length(QP.Organize.HeatVented)
    if QP.Organize.HeatVented(i)>0
        excessHeat(1,i) = x(QP.Organize.HeatVented(i));
    end
end
%pull out any dumped cooling
for i = 1:1:length(QP.Organize.CoolVented)
    if QP.Organize.CoolVented(i)>0
        excessCool(1,i) = x(QP.Organize.CoolVented(i));
    end
end
end%Ends function sortSolution