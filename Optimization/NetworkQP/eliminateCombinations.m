function [FeasibleDispatch,Cost,K] = eliminateCombinations(QP_0,K,n_best,netDemand,dt)
%This function identifies all feasible generator combinations that can meet
%the demand at this moment in time
%These feasible combinations are then tested, begining with the options
%that have the fewest number of generators
%Some combinations are avoided if they can be pre-emptively determined to
%be more costly
%The function returns the feasible generator dispatches, the cost of that 
%dispatch, and the on/off binary matrix that represents those combinations
if license('test','Distrib_Computing_Toolbox') 
    parallel = true;
else parallel = false;
end

nG = length(QP_0.constCost);
[~,n] = size(QP_0.organize);
nB = length(QP_0.Organize.Building.r);
nL = n-nG-nB;
[lines,~] = size(K);
FeasibleDispatch = zeros(lines,nG+nL+nB);
Cost = zeros(lines,1);
feasible = false(lines,1);
flag1 = zeros(lines,1);

if parallel
    parfor i = 1:lines
        [FeasibleDispatch(i,:),Cost(i), flag1(i)] = callQPsolver(QP_0,[],K(i,:));
    end
else
    nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
    [~, line] = sort(nzK); %sort the rows by number of generators that are on
    K = K(line,:);
    i = 1;
    while i<=length(K(:,1))
        [FeasibleDispatch(i,:),Cost(i), flag1(i)] = callQPsolver(QP_0,[],K(i,:));
        K = reduceK(K,netDemand,i,min(Cost(1:i)),dt);
        i = i+1;
    end
end

feasible(flag1==1) = true;
Cost = Cost+sum(ones(lines,1)*QP_0.constCost.*(K>0),2);
Cost(feasible==false) = inf;
%Cost = Cost(feasible);
%if isempty(Cost)
if nnz(feasible==true)==0
    disp('Zero feasible outcomes in eliminateCombinations: ERROR')
end
[Cost,I] = sort(Cost);
%n_best = min(n_best,length(Cost));
n_best = min(n_best,nnz(feasible==true));
K = K(I(1:n_best),:); %the n best combinations
FeasibleDispatch = FeasibleDispatch(I(1:n_best),:);

% if isfield(QP.Organize,'HeatVented')
%     HeatVent = sum(dispatch(nonzeros(QP.Organize.HeatVented),:),1)';
% else HeatVent = zeros(lines,1);%net heat loss
% end
end %ends function EliminateCombinations

function K = reduceK(K,netDemand,i,bestCost,dt)
%%reduce size of K if you don't have parallel computing toolbox
if isempty(bestCost)
    bestCost = inf;
end
if Cost < bestCost
    %remove combinations that are the best combination and additional more expensive generators
    Outs = fieldnames(netDemand);
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
            bestCostPerkWh = Cost/(netDemand.(Outs{s})*dt);
            nRow = length(K(:,1))-i;%how many rows do you have left
            include = false(nG,1);
            for j = 1:1:nG
                if QP.Organize.Dispatchable(j)
                    states = QP.Organize.States{j};
                    if any(QP.Aeq(req,states))
                        include(j) = true;
                    end
                end
            end
            if ~isempty(include) && i<length(K(:,1))
                posCheaperGen = include(logical((minRate(include)<bestCostPerkWh).*(1-(K(i,include)>0))));%possibly cheaper generators that are not on for this case, but could be on
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
            for m = 1:1:length(include)
                if K(i,include(m))==1 && i<length(K(:,1))
                    cheapGen = include(m);
                    expensiveGens = minrate(include)>minrate(cheapGen);
                    for q = 1:1:length(expensiveGens)
                        addGen = zeros(1,nG);
                        addGen(expensiveGens(q)) = 1;
                        addGen(cheapGen) = -1;
                        swapK = K(i,:)+addGen;
                        K = K~=swapK;
                    end
                end
            end
        end
    end
end
end %end function reduceK