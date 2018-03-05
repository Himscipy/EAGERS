function marginal = updateMarginalCost(Dispatch,scaleCost,dt,vH)
%%need to calculate the marginal cost of energy for any networks that have storage
global Plant
nG = length(Plant.Generator);
S = {};
marginal = [];
Type = cell(nG,1);
for i = 1:1:nG
    Type{i} = Plant.Generator(i).Type;
    if isfield(Plant.Generator(i).QPform,'Stor')
        if strcmp(Plant.Generator(i).Type,'Hydro Storage')
            S(end+1) = {'W'};
        else
            S(end+1) = fieldnames(Plant.Generator(i).QPform.output);
        end
    elseif strcmp(Plant.Generator(i).Type,'Chiller')
        if strcmp(Plant.Generator(i).Source,'Electricity')
            S(end+1) = {'E'};
        else
            S(end+1) = {'H'};
        end
    end
end

CHP = [];
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'CHP Generator')
        CHP(end+1) = i;
    end
end

if any(strcmp('E',S)) || any(strcmp('DC',S)) || any(strcmp('C',S)) || any(strcmp('Hy',S))
    include = {'CHP Generator';'Electric Generator';'Hydrogen Generator';};
    scale = scaleCost;
    if vH
        scale(:,CHP) = scale(:,CHP)*.6;
    end
    if isempty(Dispatch)
        marginal.E = mCostLoop(include,scale,[],'E',[]);
    elseif length(Dispatch(:,1))==1
        marginal.E = mCostLoop(include,scale,Dispatch,'E',[]);
    else
        marginal.E = mCostMinMax(include,scale,Dispatch,'E',[],dt);
    end
end
if any(strcmp('DC',S))
    marginal.DC = marginal.E;
end
if any(strcmp('H',S)) || any(strcmp('C',S))
    include = {'CHP Generator';'Heater';};
    scale = scaleCost;
    if vH
        scale(CHP) = scale(CHP)*.4;
    end
    if isempty(Dispatch)
        marginal.H = mCostLoop(include,scale,[],'H',[]);
    elseif length(Dispatch(:,1))==1
        marginal.H = mCostLoop(include,scale,Dispatch,'H',[]);
    else
        marginal.H = mCostMinMax(include,scale,Dispatch,'H',[],dt);
    end
end
if any(strcmp('C',S))%%chillers are unique because the cost in QP.f is zero
    include = {'Chiller'};
    if isempty(Dispatch)
        marginal.C = mCostLoop(include,scaleCost,[],'C',marginal);
    elseif length(Dispatch(:,1))==1
        marginal.C = mCostLoop(include,scaleCost,Dispatch,'C',marginal);
    else
        marginal.C = mCostMinMax(include,scaleCost,Dispatch,'C',marginal,dt);
    end
end
if any(strcmp('W',S))
    if isempty(Dispatch) || length(Dispatch(:,1))==1
        marginal.W = 1;
    else
        marginal.W.Min = 0;
        marginal.W.Max = 1;
    end
end
if any(strcmp('Hy',S)) 
    marginal.Hy = marginal.E;
%     if isempty(Dispatch) || length(Dispatch(:,1))==1
%         marginal.Hy = 1;
%     else
%         marginal.Hy.Min = 0;
%         marginal.Hy.Max = 1;
%     end
end
end%Ends function updateMarginalCost

function mc = mCostLoop(include,scale,Dispatch,S,marginal)
global Plant
nG = length(Plant.Generator);
MarginCost = [];
On = [];
for i = 1:1:nG
    if Plant.Generator(i).Enabled && any(strcmp(Plant.Generator(i).Type,include))
        if Plant.Generator(i).Enabled && strcmp(Plant.Generator(i).Type,'Chiller')
            Cratio = Plant.Generator(i).Output.Cooling(end);
            if strcmp(Plant.Generator(i).Source,'Electricity') 
                MarginCost(end+1) = marginal.E/Cratio;
            else
                MarginCost(end+1) = marginal.H/Cratio;
            end
        else
            if isempty(Dispatch)
                MarginCost(end+1) = mCost(Plant.Generator(i).QPform,[],scale(i));
            else
                MarginCost(end+1) = mCost(Plant.Generator(i).QPform,Dispatch(1,i),scale(i));
            end
        end
        if ~isempty(Dispatch) && Dispatch(1,i)>0
            On(end+1) = 1;
        else
            On(end+1) = 0;
        end
    end
    if strcmp(Plant.Generator(i).Type,'Utility') && isfield(Plant.Generator(i).QPform.output,S)
        MarginCost(end+1) = mCost(Plant.Generator(i).QPform,[],scale(i));
        On(end+1) = 1;
    end
end
if isempty(Dispatch)
    mc = mean(MarginCost);
else
    mc = max(MarginCost.*On);
    if mc == 0
        mc = min(nonzeros(MarginCost));
    end
end
end%Ends function mCostLoop

function mc = mCostMinMax(include,scale,Dispatch,S,marginal,dt)
global Plant
nG = length(Plant.Generator);
nS = length(dt);
MinC = zeros(nS,0);
MaxC = zeros(nS,0);
On = false(nS,0);
k = 0;
ChargeIndex = false(nS,1);
gen = zeros(2,0);
for i = 1:1:nG
    if Plant.Generator(i).Enabled && any(strcmp(Plant.Generator(i).Type,include))
        k = k+1;
        gen(:,end+1) = [i;k;];
        if Plant.Generator(i).Enabled && strcmp(Plant.Generator(i).Type,'Chiller')
            Cratio = Plant.Generator(i).Output.Cooling(end);
            if strcmp(Plant.Generator(i).Source,'Electricity') 
                MinC(:,k) = marginal.E.Min/Cratio;
                MaxC(:,k) = marginal.E.Max/Cratio;
            else
                MinC(:,k) = marginal.H.Min/Cratio;
                MaxC(:,k) = marginal.H.Max/Cratio;
            end
        else
            states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,1))),1);
            min_c = Plant.Generator(i).QPform.(states{1}).f(end);
            max_c = Plant.Generator(i).QPform.(states{end}).f(end) + Plant.Generator(i).QPform.(states{end}).ub(end)*Plant.Generator(i).QPform.(states{end}).H(end);
            if min_c<0
                j = 1;
                while j<length(states) && min_c<0
                    j = j+1;
                    min_c = Plant.Generator(i).QPform.(states{j}).f(end);
                end
                if min_c<0
                    min_c = 1e-3;
                end
            end
            MinC(:,k) = min_c*scale(:,i);
            MaxC(:,k) = max_c*scale(:,i);
        end
        On(:,k) = Dispatch(2:end,k)>0;
    end
    if strcmp(Plant.Generator(i).Type,'Utility') && isfield(Plant.Generator(i).QPform.output,S)
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,1))),1);
        if ~isempty(states)
            k = k+1;
            if length(states) == 1 
                MinC(:,k) = Plant.Generator(i).QPform.(states{1}).f*scale(:,i);
            else
                MinC(:,k) = Plant.Generator(i).QPform.(states{2}).f*scale(:,i);
            end
            MaxC(:,k) = Plant.Generator(i).QPform.(states{1}).f*scale(:,i);
            On(:,k) = true;
        end
    end
    if isfield(Plant.Generator(i).QPform,'Stor') && isfield(Plant.Generator(i).QPform.output,S)
        ChargeIndex = max(ChargeIndex,(Dispatch(2:end,i)-Dispatch(1:end-1,i))>0);
    end
end
if ~any(On)
    mc.Min = 0;
    mc.Max = 1;
else
    MinC(~On) = nan;
    MaxC(~On) = nan;
    GenOnDuringCharge = false(nS,k);
    ChargeIndex2 = nonzeros((1:nS)'.*ChargeIndex)+1;
    GenOnDuringCharge(ChargeIndex,gen(2,:)) = Dispatch(ChargeIndex2,gen(1,:))>0;
    if any(ChargeIndex>0) && any(any(GenOnDuringCharge))
        MaxC(GenOnDuringCharge==1) = 1.5*MaxC(GenOnDuringCharge==1); %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
    end
    minOn_t = min(MinC,[],2); 
    maxOn_t = max(MaxC,[],2);
    timedivideMin = max(1,sum(dt(minOn_t~=0)));
    timedivideMax = sum(dt(maxOn_t>0));
    minOn_t(isnan(minOn_t))=0;
    maxOn_t(isnan(maxOn_t))=0;
    mc.Min = sum(minOn_t(minOn_t~=0).*dt(minOn_t~=0))/timedivideMin;%make the cost proportional to amount of time at that cost
    mc.Max = sum(maxOn_t(maxOn_t>0).*dt(maxOn_t>0))/timedivideMax;
end
end%Ends function mCostMinMax

function mC = mCost(gen,set,scale)
states = gen.states(1:nnz(~cellfun('isempty',gen.states(:,end))),end);
mC = gen.(states{1}).f(end);
if ~isempty(set) 
    j = 1;
    while j<length(states) && gen.(states{j}).ub(end)<=set
        set = set - gen.(states{j}).ub(end);
        j = j+1;
    end
    mC = gen.(states{j}).f(end) + set*gen.(states{j}).H(end);
end
if mC<0
    j = 1;
    while j<length(states) && mC<0
        j = j+1;
        mC = gen.(states{j}).f(end);
    end
    if mC<0
        mC = 1e-3;
    end
end
mC = mC*scale;
end%Ends function mCost