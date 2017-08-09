function [GenDisp, Feasible] = FindFeasible(QPmain,Locked)
%remove ramping constraints and add them back in until it fails to isolate
%what constraint is causing it to be infeasible
%%currently unsure of what the fix is.

noRamp = QPmain.Organize.Dispatchable;
%first check it is feasible without ramping
QP = removeRamping(QPmain,noRamp);
[GenDisp, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B

% %now try fixing problem
% for i = 1:1:nnz(noRamp)
%     QP = removeRamping(QPmain,noRamp);
%     [GenDisp, Cost, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B
% end

function QP = removeRamping(QP,noRamp)
nG = length(QP.Organize.Dispatchable);
nS = length(QP.organize(:,1))-1;
ramplimit = [];
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1 && noRamp(i) ==1 %remove ramping from this generator
        ramplimit(end+1:end+2*nS,1) = (QP.Organize.Ramping(i):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.Ramping(i))';
    end
end
QP.b(ramplimit) = inf;
