function [Solution, Feasible] = FindFeasible(QP_0,Locked)
global Plant
%Try removing absorption chiller
nG = length(Plant.Generator);
AbChill = false(1,nG);
Solution.Dispatch = zeros(size(Locked));
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(Plant.Generator(i).Source,'Heat')
        AbChill(i) = true;
        Locked(:,i) = false;
    end
end
if any(AbChill)
    Locked(:,AbChill) = false;
    QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
    [x,Feasible] = callQPsolver(QP);
    if Feasible == 1
        Solution = sortSolution(x,QP);
    end
else Feasible = 0;
end
if Feasible ~=1
    %remove ramping constraints and add them back in until it fails to isolate
    %what constraint is causing it to be infeasible
    %%currently unsure of what the fix is.

    noRamp = QP_0.Organize.Dispatchable;
    %first check it is feasible without ramping
    QP_1 = removeRamping(QP_0,noRamp);
    QP = disableGenerators(QP_1,Locked,[]);%Disable generators here
    [x,Feasible] = callQPsolver(QP);
    if Feasible == 1
        Solution = sortSolution(x,QP);
    end
end
% %now try fixing problem
% for i = 1:1:nnz(noRamp)
%     QP_1 = removeRamping(QP_0,noRamp);
%     QP = disableGenerators(QP_1,Locked,[]);%Disable generators here
%     [x,Feasible] = callQPsolver(QP;
%     if Feasible == 1
%         Solution = sortSolution(x,QP);
%     end
% end
end %ends function FindFeasible

function QP = removeRamping(QP,noRamp)
nG = length(QP.Organize.Dispatchable);
nS = length(QP.organize(:,1))-1;
ramplimit = [];
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1 && noRamp(i) ==1 %remove ramping from this generator
        ramplimit(end+1:end+nS,1) = (QP.Organize.Ramping(i):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.Ramping(i))';
    end
end
QP.b(ramplimit) = inf;
QP.b(ramplimit+1) = inf;
end%Endss function removeRamping