function Solution = sortSolution(x,QP)
[m,n] = size(QP.organize);
nS = m-1;
nG = length(QP.constCost);
nB = length(QP.Organize.Building.r);
nL = n-nG-nB;
nH = nnz(QP.Organize.Hydro);
Solution.Dispatch = zeros(nS+1,nG);
Solution.hydroSOC = zeros(nS,nH);
Solution.excessHeat = zeros(nS,nnz(QP.Organize.HeatVented));
Solution.LineFlows = zeros(nS,nL);
Solution.LineLoss = zeros(nS,nL);
Solution.Buildings = zeros(nS,nB);
for i = 1:1:nG
    if isfield(QP,'Renewable') && any(QP.Renewable(:,i)~=0)
        Solution.Dispatch(2:end,i) = QP.Renewable(:,i);
    else
        Out_vs_State = QP.Organize.Out_vs_State{1,i};%linear maping between state and output
        for t = 1:1:nS+1
            if ~isempty(QP.organize{t,i})
                states = QP.organize{t,i};
                for j = 1:1:length(states)
                    Solution.Dispatch(t,i) = Solution.Dispatch(t,i) + Out_vs_State(j)*x(states(j)); %record this combination of outputs (SOC for storage)
                end
            end
        end
    end
    if QP.Organize.Hydro(i)
        for t = 1:1:nS
            %Get SOC of each generator into a matrix for all time steps
            Solution.hydroSOC(t,i) = x(QP.organize{t+1,i}+1,1);
        end
    end
end
for i = 1:1:nL
    for t = 1:1:nS
        Solution.LineFlows(t,i) = sum(x(QP.organize{t+1,i+nG}));%power transfer
        if length(QP.Organize.States{i+nG})>1
            Solution.LineLoss(t,i) = sum(x(QP.organize{t+1,i+nG}+1)); %down (positive) lines
            Solution.LineLoss(t,i) = Solution.LineLoss(t,i) + sum(x(QP.organize{t+1,i+nG}+2)); %up (negative) lines
        end
    end
end
for i = 1:1:nB
    for t = 1:1:nS
        Solution.Buildings(t,i) = x(QP.organize{t+1,i+nG+nL},1);
    end
end

Solution.Dispatch(abs(Solution.Dispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
%pull out any dumped heat
for i = 1:1:length(QP.Organize.HeatVented)
    if QP.Organize.HeatVented(i)>0
        for t = 1:1:nS
            Solution.excessHeat(t,i) = x(QP.Organize.HeatVented(i)+ QP.Organize.t1States*(t-1));
        end
    end
end
end%Ends function sortSolution