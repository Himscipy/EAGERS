function marginal = updateMarginalCost(Dispatch,scaleCost,dt,fit)
global Plant
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
networkNames = networkNames(~strcmp('Location',networkNames));
for i = 1:1:length(networkNames)
    if strcmp(networkNames{i},'Electrical')
        Out.E = [];
    elseif strcmp(networkNames{i},'DistrictHeat')
        Out.H = [];
    elseif strcmp(networkNames{i},'DistrictCool')
        Out.C = [];
    elseif strcmp(networkNames{i},'Hydro')
        Out.W = [];
    end
end
storType = Out;
marginal = Out;
nG = length(Plant.Generator);
nS = length(Dispatch(:,1))-1;
UB = zeros(nG,1);
marginCost = zeros(2,nG);
CHP = [];
I = false(1,nG);
for i = 1:1:nG
    if Plant.Generator(i).Enabled && ~isempty(Plant.Generator(i).QPform.states)
        
        if isfield(Plant.Generator(i).QPform,'Stor')
            if strcmp(Plant.Generator(i).Source,'Electricity')
                storType.E(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                storType.H(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                storType.C(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Water')
                storType.W(end+1) = i;
            end
        else
            I(i) = true;
            if strcmp(Plant.Generator(i).Type,'Utility')
                states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,1))),1);
                if ~isempty(states)
                    for j = 1:1:length(states)
                        UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(1);
                    end
                    if length(states) == 1 
                        marginCost(1,i) = Plant.Generator(i).QPform.(states{1}).f;
                    else marginCost(1,i) = Plant.Generator(i).QPform.(states{2}).f;
                    end
                    marginCost(2,i) = Plant.Generator(i).QPform.(states{1}).f;
                end
            else                
                states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,fit))),fit);
                if~isempty(states)
                    for j = 1:1:length(states)
                        UB(i) = UB(i) + Plant.Generator(i).QPform.(states{j}).ub(fit);
                    end
                    marginCost(1,i) = Plant.Generator(i).QPform.(states{1}).f(fit);
                    marginCost(2,i) = Plant.Generator(i).QPform.(states{end}).f(fit) + Plant.Generator(i).QPform.(states{end}).ub(fit)*Plant.Generator(i).QPform.(states{end}).H(fit);
                    if marginCost(1,i)<0
                        j = 1;
                        while j<length(states) && marginCost(1,i)<0
                            j = j+1;
                            marginCost(1,i) = Plant.Generator(i).QPform.(states{j}).f(fit);
                        end
                        if marginCost(1,i)<0
                            marginCost(1,i) = 1e-3;
                        end
                    end
                end
            end
        end
        if I(i)
            S = fieldnames(Plant.Generator(i).QPform.output);
            if length(S)==2 && ismember('H',S) && ismember('E',S)
                Out.E(end+1) = i;
                CHP(end+1) = i;
            else
                Out.(S{1})(end+1) = i;
            end
        end
    end
end
GenOnDuringCharge = false(nS,nG);
minMarginCost = nan(nS,nG);
maxMarginCost = nan(nS,nG);

for i = 1:1:nG
    if I(i) && UB(i)>0
        minMarginCost(:,i) = scaleCost(:,i).*marginCost(1,i);
        maxMarginCost(:,i) = scaleCost(:,i).*marginCost(2,i);
    end
end
Type = fieldnames(storType);
for i = 1:1:length(Type)
    gen = Out.(Type{i});
    if strcmp(Type{i},'C') %chiller cost shows up as electric or heating load
        Capacity = 0;
        DispatchedC = 0;
        abChill = [];
        for j = 1:1:length(gen)
            Cratio = Plant.Generator(gen(j)).Output.Cooling(end);
            Capacity = Capacity + nS*Plant.Generator(gen(j)).Size;
            DispatchedC = DispatchedC + sum(Dispatch(:,gen(j)));
            if strcmp(Plant.Generator(gen(j)).Source,'Electricity') 
                minMarginCost(:,gen(j)) = max(0,min(minMarginCost(:,Out.E),[],2))/Cratio;
                maxMarginCost(:,gen(j)) = max(maxMarginCost(:,Out.E),[],2)/Cratio;
            else
                abChill(end+1) = gen(j);
                minMarginCost(:,gen(j)) = min(minMarginCost(:,Out.H),[],2)/Cratio;
                maxMarginCost(:,gen(j)) = max(maxMarginCost(:,Out.H),[],2)/Cratio;
            end
        end
        %%Need to replace with calculation for excessheat
        if ~isempty(abChill) && DispatchedC<0.2*Capacity
            minMarginCost(:,abChill) = -0.5*max(maxMarginCost(:,gen),[],2);%0;
            maxMarginCost(:,gen) = 0.5*maxMarginCost(:,gen);
        end
    end
    
    ChargeIndex = false(nS,1);
    for j = 1:1:length(storType.(Type{i}))
        ChargeIndex = max(ChargeIndex,Dispatch(2:end,storType.(Type{i})(j))-Dispatch(1:end-1,storType.(Type{i})(j))>0);%[false;(Dispatch(2:end,stor(i))-Dispatch(1:end-1,stor(i)))>0]; %timesteps where charging
    end
    GenOnDuringCharge(ChargeIndex,gen) = Dispatch(ChargeIndex,gen)>0;
    if any(ChargeIndex>0) && any(any(GenOnDuringCharge))
        maxMarginCost(GenOnDuringCharge==1) = 1.5*maxMarginCost(GenOnDuringCharge==1); %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
    end
    MinThisType = minMarginCost(:,gen);
    MaxThisType = maxMarginCost(:,gen);
    if ~isempty(CHP)
        if strcmp(Type{i},'H')
            MinThisType(:,end+1:end+length(CHP)) = minMarginCost(:,CHP)*.4; %assign 40% of the generator cost to the heat production
            MaxThisType(:,end+1:end+length(CHP)) = maxMarginCost(:,CHP)*.4; %assign 40% of the generator cost to the heat production
        end
    end
    if isempty(MinThisType)
        %%this is creating a marginal cost for hydro storage only dispatch
        marginal.(Type{i}).Min = 0;
        marginal.(Type{i}).Max = 1;
    else
        minOn_t = min(MinThisType,[],2); 
        maxOn_t = max(MaxThisType,[],2);
        timedivideMin = max(1,sum(dt(minOn_t~=0)));
        timedivideMax = sum(dt(maxOn_t>0));
        minOn_t(isnan(minOn_t))=0;
        maxOn_t(isnan(maxOn_t))=0;
        marginal.(Type{i}).Min = sum(minOn_t(minOn_t~=0).*dt(minOn_t~=0))/timedivideMin;%make the cost proportional to amount of time at that cost
        marginal.(Type{i}).Max = sum(maxOn_t(maxOn_t>0).*dt(maxOn_t>0))/timedivideMax;
    end
end