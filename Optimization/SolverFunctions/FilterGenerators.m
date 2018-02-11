function Solution = FilterGenerators(QP_0,Solution,Locked,Timestamp)
global Plant OnOff
dt = (Timestamp(2:end) - Timestamp(1:end-1))*24;
nS = length(dt);
dGen = find(QP_0.Organize.Dispatchable==1);
%% load parameters
nG = length(Plant.Generator);
dX = zeros(nS,nG);
UB = zeros(nG,1);
for i = 1:1:nG
    if ~isempty(Plant.Generator(i).QPform.states)
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(end);
        end
    end
    dX(:,i) = dt*Plant.Generator(i).VariableStruct.dX_dt;
end
%% current best guess
FirstDisp = Solution.Dispatch;
Cost = sum(NetCostCalc(Solution.Dispatch,Timestamp,'Dispatch'));
%% Proceed with rules
index = (1:nS)';
for j = 1:1:length(dGen)
    %% Find when starts and stops occur in optimal dispatch
    i = dGen(j);
    starts = nonzeros(index.*(((FirstDisp(2:end,i)>0)-(FirstDisp(1:nS,i)>0))>0)); % if on at t = 3, start = 3, where IC is t=0
    if ~OnOff(i)&& FirstDisp(1,i)>0
        starts = [1;starts];
    end
    stops = nonzeros(index.*(((FirstDisp(1:nS,i)>0)-(FirstDisp(2:end,i)>0))>0)); % if off at t = 12, stop = 12, where IC is t=0
    if ~Locked(1,i) || ~OnOff(i)
        %% Rule 1: turn off for longer time at start if possible
        if~isempty(starts)
            p = 0;
            RampUp = dX(starts(1),i);
            while RampUp<FirstDisp(starts(1)+1,i) && (starts(1)-p>0)
                RampUp = RampUp+dX(starts(1)-p,i);
                p = p+1;
            end
            if starts(1)-p>0
                L2 = Locked;
                L2(1:(starts(1)-p+1),i) = false;
                QP = disableGenerators(QP_0,L2,[]);%Disable generators here
                [x,Feasible] = callQPsolver(QP);
                if Feasible == 1
                    SolNew = sortSolution(x,QP);
                    newCost = sum(NetCostCalc(SolNew.Dispatch,Timestamp,'Dispatch'));
                    if newCost<Cost
                        Locked = L2;
                        Cost = newCost;
                        Solution = SolNew;
                    end
                end
            end
            if length(starts)>1
                starts = starts(2:end);
            else starts =[];
            end
        end
    end
    %% Rule 2: If off for a long enough segment in optimal dispatch, try turning off for as much of that as possible given ramp rates
    for k = 1:1:length(starts)
        if ~isempty(stops) && length(stops)>=k
            if sum(dX(stops(k):starts(k)-1,i))>(FirstDisp(stops(k),i) + FirstDisp(starts(k)+1,i)) %can ramp all the way down, and back up, and be off for 1 step
                L2 = Locked;
                %find step when it can hit zero given setting at Disp(stops(k)
                n=1;
                RampDown = dX(stops(k),i);
                while RampDown<FirstDisp(stops(k),i)
                    RampDown = RampDown+dX(stops(k)+n,i);
                    n = n+1;
                end
                p = 1;
                RampUp = dX(starts(k),i);
                while RampUp<FirstDisp(starts(k)+1,i)
                    RampUp = RampUp+dX(starts(k)-p,i);
                    p = p+1;
                end
                L2((stops(k)+n):(starts(k)-p+1),i) = false;
                QP = disableGenerators(QP_0,L2,[]);%Disable generators here
                [x,Feasible] = callQPsolver(QP);
                if Feasible == 1
                    SolNew = sortSolution(x,QP);
                    newCost = sum(NetCostCalc(SolNew.Dispatch,Timestamp,'Dispatch'));
                    if newCost<Cost
                        Locked = L2;
                        Cost = newCost;
                        Solution = SolNew;
                    end
                end
            end
        end
    end
    %% Rule 3: try turning off @ end
    if length(stops)>length(starts)
        n=stops(end);
        RampDown = dX(n,i);
        while RampDown<FirstDisp(stops(end),i) && n<(nS) && n>0
            RampDown = RampDown+dX(n,i);
            n = n-1;
        end
        if n<(nS)
            L2 = Locked;
            L2((n+1):nS+1,i) = false;
            QP = disableGenerators(QP_0,L2,[]);%Disable generators here
            [x,Feasible] = callQPsolver(QP);
            if Feasible == 1
                SolNew = sortSolution(x,QP);
                newCost = sum(NetCostCalc(SolNew.Dispatch,Timestamp,'Dispatch'));
                if newCost<Cost
                    Locked = L2;
                    Cost = newCost;
                    Solution = SolNew;
                end
            end
        end
    end
end

for j = 1:1:length(dGen)
    %% Find when starts and stops occur in current Locked matrix
    i = dGen(j);
    starts = nonzeros(index.*((Locked(2:end,i)-Locked(1:nS,i))>0)); % if on at t = 3, start = 3, where IC is t=0
    stops = nonzeros(index.*((Locked(1:nS,i)-Locked(2:end,i))>0)); % if off at t = 12, stop = 12, where IC is t=0

    %% Rule 4: if on for a short time and sufficient capacity in remaining active generators, turn gen off completely
    if Locked(1,i)
        if length(stops)>1
            stops = stops(2:end);
        else stops =[];
        end
    end
    for k = 1:1:length(stops)
        if sum(dX(starts(k):stops(k),i))<UB(i) && sum(Locked(:,i))<floor(nS/4)%can only ramp to 1/2 power and less than 1/4 of horizon
            L2 = Locked;
            L2(starts(k):(stops(k)+1),i)= false;
            QP = disableGenerators(QP_0,L2,[]);%Disable generators here
            [x,Feasible] = callQPsolver(QP);
            if Feasible == 1
                SolNew = sortSolution(x,QP);
                newCost = sum(NetCostCalc(SolNew.Dispatch,Timestamp,'Dispatch'));
                if newCost<Cost
                    Locked = L2;
                    Cost = newCost;
                    Solution = SolNew;
                end
            end
        end
    end
end
% % prevent flicker
% for j = 1:1:length(dGen)
%     i = dGen(j);
%     starts = nonzeros(index.*((Locked(2:end,i)-Locked(1:nS,i))>0));
%     stops = nonzeros(index.*((Locked(1:nS,i)-Locked(2:end,i))>0));
%     if ~isempty(starts) && ~isempty(stops)
%         stops = stops(stops<starts(end));
%     end
%     if ~isempty(stops) && ~isempty(stops)
%         starts = starts(starts>stops(1));
%         for k = 1:1:length(starts)
%             if starts(k)==stops(k)+1 %if you are starting immediately after stopping
%                 L2 = Locked;
%                 L2(stops(k)+1,i) = true;
%                 QP = disableGenerators(QP_0,L2,[]);%Disable generators here
%                 [x,Feasible] = callQPsolver(QP);
%                 if Feasible == 1
%                     SolNew = sortSolution(x,QP);
%                     newCost = sum(NetCostCalc(SolNew.Dispatch,Timestamp,'Dispatch'));
%                     if newCost<Cost
%                         Locked = L2;
%                         Cost = newCost;
%                         Solution = SolNew;
%                     end
%                 end
%             end
%         end
%     end
% end