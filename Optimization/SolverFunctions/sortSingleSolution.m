function [Dispatch,LineLoss,excessHeat,hydroSOC] = sortSingleSolution(x,QP)
[~,n] = size(QP.organize);
nG = length(QP.constCost);
nB = length(QP.Organize.Building.r);
nL = n-nG-nB;
nH = nnz(QP.Organize.Hydro);
hydroSOC = zeros(1,nH);
excessHeat = zeros(1,nnz(QP.Organize.HeatVented));
LineFlows = zeros(1,nL);
LineLoss = zeros(1,nL);
Buildings = zeros(1,nB);
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
    if QP.Organize.Hydro(i) == 1
        hydroSOC(1,i) = x(QP.organize{1,i}+1,1);
    end
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
    Buildings(1,i) = x(QP.organize{1,i+nG+nL},1);
end
Dispatch = [GenDisp,LineFlows,Buildings];
%pull out any dumped heat
for i = 1:1:length(QP.Organize.HeatVented)
    if QP.Organize.HeatVented(i)>0
        excessHeat(1,i) = x(QP.Organize.HeatVented(i));
    end
end
end%Ends function sortSolution