function Solution = cQP_Feasibility(OptimalState,Forecast,scaleCost,Date)
global Plant
dt = (Date(2:end) - Date(1:end-1))*24;
marginCost = updateMarginalCost(OptimalState,scaleCost,dt,2);
QP_0 = updateMatrices(Plant.OpMatB,Date,scaleCost,marginCost,Forecast,[]); %update fit B matrices
nS = length(dt);
nG = length(Plant.Generator);
Locked = true(nS+1,nG);
PossiblyOn = false(nS,nG);
for i = 1:1:nG
    if ~Plant.Generator(i).Enabled
        Locked(:,i) = 0;
    end
end
LB = zeros(1,nG);
for i = 1:1:nG
    if QP_0.Organize.Dispatchable(i) ==1
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            LB(i) = LB(i) + Plant.Generator(i).QPform.(states{j}).lb(2);
        end
        Locked(OptimalState(:,i)<LB(i),i)=false;%Default to off when initial dispatch is below LB
        PossiblyOn(OptimalState(2:end,i)>0 & OptimalState(2:end,i)<LB(i)) = true;
    end
end
%% make sure it can shut down in time from initial condition
for i = 1:1:nG
    if QP_0.Organize.Dispatchable(i) ==1
        if OptimalState(1,i)>0 && ~all(Locked(:,i))
            r = QP_0.Organize.Ramping(i)+1;
            D = OptimalState(1,i);
            t = 1;
            while D>0
                D = D - QP_0.b(r);
                if D>0 && ~Locked(t+1,i)
                    Locked(t+1,i) = true;
                end
                t = t+1;
                r = r+QP_0.Organize.t1ineq;
            end
        end
    end
end

Feasible = 0;
attempt = 0;
while Feasible ~= 1 && attempt <1
    [x,Feasible] = checkFeas(QP_0,Locked);
    %% add logic here adjusting Locked to check other conditions if it comes up infeasible
    %PossiblyOn is a 24 x nG matrix of when things could be on, note Locked is 25 x nG
    attempt = 1;
end
%%%%

if Feasible==1
    QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
    Solution = sortSolution(x,QP);
    Solution = FilterGenerators(QP_0,Solution,Locked,Date);
    %add something to see if FilterGenerators changed anything
else
    disp('error: Cannot Find Feasible Dispatch');
end
end%ends function cQP_Feasibility

function [x,Feasible] = checkFeas(QP_0,Locked)
QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
[x,Feasible] = callQPsolver(QP);%this is the dispatch with fit B
end%Ends function checkFeas