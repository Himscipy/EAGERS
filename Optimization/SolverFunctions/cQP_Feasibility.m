function [Solution,LBrelax] = cQP_Feasibility(OptimalState,Forecast,scaleCost,Date)
global Plant
dt = (Date(2:end) - Date(1:end-1))*24;
marginCost = updateMarginalCost(OptimalState,scaleCost,dt,[]);
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
        %PossiblyOn is a 24 x nG matrix of when things could be on, note Locked is 25 x nG
        Locked(OptimalState(:,i)<=LB(i),i) = false;%Default to on unless offline in initial dispatch
        %Locked(OptimalState(:,i)<LB(i),i)=false;%Default to off when initial dispatch is below LB
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
LBrelax = 1;

while Feasible ~= 1 && attempt <6
    %attempt: integer value describing the number of attempts before reaching
    %feasibility, this determines how close components must be to their lower
    %bound from below to be considered online
    %n represents the percent of lower bounds 
    %on your first try, just use the locked matrix given, then do unit
    %commitment based on OptimalState>LB*percLB
    percLB = [0.9, 0.75, 0.5, 0.2, 0.1, 0, -1];

    if attempt>0%second try, lower limit for online threshold
        LBrelax = percLB(attempt);
        %only change label for unit commitment gens, and don't change the
        %label for initial conditions
        for i = 1:1:nG
            if QP_0.Organize.Dispatchable(i) ==1
                Locked(2:end,i) = (OptimalState(2:end,i)>(LB(i)*LBrelax));%Default to on unless offline in initial dispatch
            end
        end
    end
    [x,Feasible] = checkFeas(QP_0,Locked);
    attempt = attempt+1;
end


if Feasible==1
    QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
    Solution = sortSolution(x,QP);
    Solution = FilterGenerators(QP_0,Solution,Locked,Date);
    %add something to see if FilterGenerators changed anything
else
    disp('error: Cannot Find Feasible Dispatch');
end
Solution.LBrelax = LBrelax;
end%ends function cQP_Feasibility

function [x,Feasible] = checkFeas(QP_0,Locked)
%% this function finds a feasible solution for the cQP method
% it is called by cQP_Feasibility
%% inputs: QP: structure of non-sparse matrices for quadratic programming
%Locked: boolean matrix (nS+1) x nG describing which components are online,
%the first row should align with the CurrentState.GenDisp, because the first row is
%the initial conditions, non-unit commitment components should always be
%online (true, 1 classification)

%% outputs: x: real valued vector of setpoints for optimal dispatch, this vector 
%is sorted into a dispatch solution in cQP_Feasibility
%Feasible: this is 1 if the solution is feasible

QP = disableGenerators(QP_0,Locked,[]);%Disable generators here
[x,Feasible] = callQPsolver(QP);


end%Ends function checkFeas