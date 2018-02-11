function [BestDispatch,Alt] = eliminateCombinations(QP_0,K,Alt,IC,netDemand,dt,t)
%This function identifies all feasible generator combinations that can meet
%the demand at this moment in time
%These feasible combinations are then tested, begining with the options
%that have the fewest number of generators
%Some combinations are avoided if they can be pre-emptively determined to
%be more costly
%The function returns the feasible generator dispatches, the cost of that 
%dispatch, and the on/off binary matrix that represents those combinations
[lines,~] = size(K);
if lines>1 && license('test','Distrib_Computing_Toolbox') 
    parallel = true;
else
    parallel = false;
end

nG = length(QP_0.constCost);
[~,n] = size(QP_0.organize);
nB = length(QP_0.Organize.Building.r);
nL = n-nG-nB;
nH = nnz(QP_0.Organize.Hydro);

Dispatch = zeros(lines,n);
LineLoss = zeros(lines,nL);
excessHeat = zeros(lines,nnz(QP_0.Organize.HeatVented));
excessCool = zeros(lines,nnz(QP_0.Organize.CoolVented));
hydroSOC = zeros(lines,nH);
Temperature = zeros(lines,nB);
Heating = zeros(lines,nB);
Cooling = zeros(lines,nB);
Cost = zeros(lines,1);
feasible = false(lines,1);
flag1 = zeros(lines,1);

if parallel
    parfor i = 1:lines
        QP = disableGenerators(QP_0,[],K(i,:));%Disable generators here
        [x, flag1(i)] = callQPsolver(QP);
        if flag1(i)==1
            Cost(i) = 0.5*x'*QP.H*x + x'*QP.f;
            [Dispatch(i,:),LineLoss(i,:),excessHeat(i,:),excessCool(i,:),hydroSOC(i,:),Temperature(i,:),Heating(i,:),Cooling(i,:)] = sortSingleSolution(x,QP);
        end
    end
    Cost = Cost+sum(ones(lines,1)*QP_0.constCost.*(K>0),2);
else
    nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
    [~, line] = sort(nzK); %sort the rows by number of generators that are on
    K = K(line,:);
    i = 1;
    while i<=length(line)
        QP = disableGenerators(QP_0,[],K(i,:));%Disable generators here
        [x, flag1(i)] = callQPsolver(QP);
        if flag1(i)==1
            Cost(i) = 0.5*x'*QP.H*x + x'*QP.f;
            Cost(i) = Cost(i)+sum(QP_0.constCost.*(K(i,:)>0),2);
            [Dispatch(i,:),LineLoss(i,:),excessHeat(i,:),excessCool(i,:),hydroSOC(i,:),Temperature(i,:),Heating(i,:),Cooling(i,:)] = sortSingleSolution(x,QP);
        end
        if i<length(line)
            K = reduceK(K,QP,netDemand,i,Cost(i),min(Cost(1:i-1)),dt,t);
        end
        i = i+1;
    end
end

feasible(flag1==1) = true;
Cost(feasible==false) = inf;
[Cost,I] = sort(Cost);
nP = nnz(feasible==true);
K = K(I(1:nP),:); %the n best combinations

Alt.Binary{t} = K>0;
Alt.Disp{t} = Dispatch(I(1:nP),:);
Alt.Cost{t} = Cost(1:nP)-Cost(1);
Alt.Generators{t} = Dispatch(I(1:nP),1:nG);
Alt.LineFlows{t} = Dispatch(I(1:nP),nG+1:nG+nL);
Alt.Buildings{t} = Dispatch(I(1:nP),nG+nL+1:nG+nL+nB);
Alt.LineLoss{t} = LineLoss(I(1:nP),:);
Alt.excessHeat{t} = excessHeat(I(1:nP),:);
Alt.excessCool{t} = excessCool(I(1:nP),:);
Alt.hydroSOC{t} = hydroSOC(I(1:nP),:);
Alt.Temperature{t} = Temperature(I(1:nP),:);
Alt.Heating{t} = Heating(I(1:nP),:);
Alt.Cooling{t} = Cooling(I(1:nP),:);
if isempty(Alt.Disp{t})
    disp(['No feasible combination of generators at step' num2str(t)]);
    BestDispatch = IC;
else
    BestDispatch = Alt.Disp{t}(1,:);
end
end %ends function EliminateCombinations

function K = reduceK(K,QP,netDemand,i,Cost,bestCost,dt,t)
%%reduce size of K if you don't have parallel computing toolbox
if isempty(bestCost)
    bestCost = inf;
end
if Cost < bestCost
    %remove combinations that are the best combination and additional more expensive generators
    Outs = fieldnames(netDemand);
    nG = length(QP.constCost);
    for s = 1:1:length(Outs)
        req = [];
        if strcmp(Outs{s},'E')
            req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
        elseif strcmp(Outs{s},'H')
            req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
        elseif strcmp(Outs{s},'C')
            req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
        elseif strcmp(Outs{s},'W')
            req = QP.Organize.Balance.Hydro; %rows of Aeq associated with hydro demand
        end

        if length(req) ==1 % The following definitely works with only 1 node
            bestCostPerkWh = Cost/(netDemand.(Outs{s})(t)*dt);
            nRow = length(K(:,1))-i;%how many rows do you have left
            include = false(1,nG);
            minRate = inf+zeros(1,nG);
            for j = 1:1:nG
                if QP.Organize.Dispatchable(j)
                    states = QP.Organize.States{j};
                    if any(QP.Aeq(req,states))
                        include(j) = true;
                        minRate(j) = QP.f(states(1));
                        if minRate(j)== 0 && QP.Aeq(QP.Organize.Balance.Electrical,states(1))<0
                            minRate(j) = -.1*QP.Aeq(QP.Organize.Balance.Electrical,states(1));
                        elseif minRate(j)== 0 && QP.Aeq(QP.Organize.Balance.DistrictHeat,states(1))<0
                            minRate(j) = -.1*QP.Aeq(QP.Organize.Balance.DistrictHeat,states(1));
                        end
                    end
                end
            end
            if i<length(K(:,1))
                posCheaperGen = nonzeros((1:nG).*include.*((minRate<bestCostPerkWh) & (K(i,:)==0)));%possibly cheaper generators that are not on for this case, but could be on
            else
                posCheaperGen = [];
            end
            if ~isempty(posCheaperGen)
                ckRows = ~any((ones(nRow,1)*K(i,include)-K(i+1:end,include))>0,2); %identify future rows that have the current set of active generators + more
                noCheaper = ~any(K(i+1:end,posCheaperGen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
                rmRows = ckRows&noCheaper;
                if nnz(rmRows)>0
                    K = K([true(i,1);~rmRows],:);  
                end
            end
            %remove combinations that swap the current
            %combination with a more expensive generator
            for m = 1:1:nG
                if include(m) && K(i,m)>0 && i<length(K(:,1))
                    expensiveGens = nonzeros((1:nG).*include.*(minRate>minRate(m)));
                    for q = 1:1:length(expensiveGens)
                        if K(i,expensiveGens(q))==0
                            addGen = zeros(1,nG);
                            addGen(expensiveGens(q)) = expensiveGens(q);
                            addGen(m) = -m;
                            swapK = K(i,:)+addGen;
                            keep = ~ismember(K,swapK,'rows');
                            K = K(keep,:);
                        end
                    end
                end
            end
        end
    end
end
end %end function reduceK