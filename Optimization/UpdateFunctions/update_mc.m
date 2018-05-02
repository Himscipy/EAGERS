function marginal = update_mc(gen,dispatch,scale_cost,dt,v_h)
%%need to calculate the marginal cost of energy for any networks that have storage
n_g = length(gen);
s = {};
marginal = [];
type = cell(n_g,1);
for i = 1:1:n_g
    type{i} = gen(i).Type;
    if isfield(gen(i).QPform,'Stor')
        if strcmp(gen(i).Type,'Hydro Storage')
            s(end+1) = {'W'};
        else
            s(end+1) = fieldnames(gen(i).QPform.output);
        end
    elseif strcmp(gen(i).Type,'Chiller')
        if strcmp(gen(i).Source,'Electricity')
            s(end+1) = {'E'};
        else
            s(end+1) = {'H'};
        end
    end
end

chp = [];
for i = 1:1:n_g
    if strcmp(gen(i).Type,'CHP Generator')
        chp(end+1) = i;
    end
end
if ~isempty(chp) && isempty(v_h)
    v_h = true;
end
    
if any(strcmp('E',s)) || any(strcmp('DC',s)) || any(strcmp('C',s)) || any(strcmp('Hy',s))
    include = {'CHP Generator';'Electric Generator';'Hydrogen Generator';};
    scale = scale_cost;
    if v_h
        scale(:,chp) = scale(:,chp)*.7;
    end
    if isempty(dispatch)
        marginal.E = m_loop(gen,include,scale,[],'E',[]);
    elseif length(dispatch(:,1))==1
        marginal.E = m_loop(gen,include,scale,dispatch,'E',[]);
    else
        marginal.E = min_max_cost(gen,include,scale,dispatch,'E',[],dt);
    end
end
if any(strcmp('DC',s))
    marginal.DC = marginal.E;
end
if any(strcmp('H',s)) || any(strcmp('C',s))
    include = {'CHP Generator';'Heater';};
    scale = scale_cost;
    if v_h
        scale(:,chp) = scale(:,chp)*.3;
    end
    if isempty(dispatch)
        marginal.H = m_loop(gen,include,scale,[],'H',[]);
    elseif length(dispatch(:,1))==1
        marginal.H = m_loop(gen,include,scale,dispatch,'H',[]);
    else
        marginal.H = min_max_cost(gen,include,scale,dispatch,'H',[],dt);
    end
end
if any(strcmp('C',s))%%chillers are unique because the cost in QP.f is zero
    include = {'Chiller'};
    if isempty(dispatch)
        marginal.C = m_loop(gen,include,scale_cost,[],'C',marginal);
    elseif length(dispatch(:,1))==1
        marginal.C = m_loop(gen,include,scale_cost,dispatch,'C',marginal);
    else
        marginal.C = min_max_cost(gen,include,scale_cost,dispatch,'C',marginal,dt);
    end
end
if any(strcmp('W',s))
    if isempty(dispatch) || length(dispatch(:,1))==1
        marginal.W = 1;
    else
        marginal.W.Min = 0;
        marginal.W.Max = 1;
    end
end
if any(strcmp('Hy',s)) 
    marginal.Hy = marginal.E;
%     if isempty(Dispatch) || length(Dispatch(:,1))==1
%         marginal.Hy = 1;
%     else
%         marginal.Hy.Min = 0;
%         marginal.Hy.Max = 1;
%     end
end
end%Ends function update_mc

function mc = m_loop(gen,include,scale,dispatch,s,marginal)
n_g = length(gen);
margin_cost = [];
on = [];
for i = 1:1:n_g
    if gen(i).Enabled && any(strcmp(gen(i).Type,include))
        if gen(i).Enabled && strcmp(gen(i).Type,'Chiller')
            Cratio = gen(i).Output.Cooling(end);
            if strcmp(gen(i).Source,'Electricity') 
                margin_cost(end+1) = marginal.E/Cratio;
            else
                margin_cost(end+1) = marginal.H/Cratio;
            end
        else
            if isempty(dispatch)
                margin_cost(end+1) = m_cost(gen(i).QPform,[],scale(i));
            else
                margin_cost(end+1) = m_cost(gen(i).QPform,dispatch(1,i),scale(i));
            end
        end
        if ~isempty(dispatch) && dispatch(1,i)>0
            on(end+1) = 1;
        else
            on(end+1) = 0;
        end
    end
    if strcmp(gen(i).Type,'Utility') && isfield(gen(i).QPform.output,s)
        margin_cost(end+1) = m_cost(gen(i).QPform,[],scale(i));
        on(end+1) = 1;
    end
end
if isempty(dispatch)
    mc = mean(margin_cost);
else
    mc = max(margin_cost.*on);
    if mc == 0
        mc = min(nonzeros(margin_cost));
    end
end
end%Ends function m_loop

function mc = min_max_cost(gen,include,scale,dispatch,s,marginal,dt)
n_g = length(gen);
n_s = length(dt);
min_c = zeros(n_s,0);
max_c = zeros(n_s,0);
on = false(n_s,0);
k = 0;
charge_index = false(n_s,1);
gen_index = zeros(2,0);
for i = 1:1:n_g
    if gen(i).Enabled && any(strcmp(gen(i).Type,include))
        k = k+1;
        gen_index(:,end+1) = [i;k;];
        if gen(i).Enabled && strcmp(gen(i).Type,'Chiller')
            Cratio = gen(i).Output.Cooling(end);
            if strcmp(gen(i).Source,'Electricity') 
                min_c(:,k) = marginal.E.Min/Cratio;
                max_c(:,k) = marginal.E.Max/Cratio;
            else
                min_c(:,k) = marginal.H.Min/Cratio;
                max_c(:,k) = marginal.H.Max/Cratio;
            end
        else
            states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,1))),1);
            min_c_k = gen(i).QPform.(states{1}).f(end);
            max_c_k = gen(i).QPform.(states{end}).f(end) + gen(i).QPform.(states{end}).ub(end)*gen(i).QPform.(states{end}).H(end);
            if min_c_k<0
                j = 1;
                while j<length(states) && min_c_k<0
                    j = j+1;
                    min_c_k = gen(i).QPform.(states{j}).f(end);
                end
                if min_c_k<0
                    min_c_k = 1e-3;
                end
            end
            min_c(:,k) = min_c_k*scale(:,i);
            max_c(:,k) = min(3*min_c_k,max_c_k)*scale(:,i);
        end
        on(:,k) = dispatch(2:end,i)>0;
    end
    if strcmp(gen(i).Type,'Utility') && isfield(gen(i).QPform.output,s)
        states = gen(i).QPform.states(1:nnz(~cellfun('isempty',gen(i).QPform.states(:,1))),1);
        if ~isempty(states)
            k = k+1;
            if length(states) == 1 
                min_c(:,k) = gen(i).QPform.(states{1}).f*scale(:,i);
            else
                min_c(:,k) = gen(i).QPform.(states{2}).f*scale(:,i);
            end
            max_c(:,k) = gen(i).QPform.(states{1}).f*scale(:,i);
            on(:,k) = true;
        end
    end
    if isfield(gen(i).QPform,'Stor') && isfield(gen(i).QPform.output,s)
        charge_index = max(charge_index,(dispatch(2:end,i)-dispatch(1:end-1,i))>0);
    end
end
if ~any(on)
    mc.Min = 0;
    mc.Max = 1;
else
    min_c(~on) = nan;
    max_c(~on) = nan;
    gen_on_during_charge = false(n_s,k);
    charge_index2 = nonzeros((1:n_s)'.*charge_index)+1;
    gen_on_during_charge(charge_index,gen_index(2,:)) = dispatch(charge_index2,gen_index(1,:))>0;
    if any(charge_index>0) && any(any(gen_on_during_charge))
        max_c(gen_on_during_charge==1) = 1.5*max_c(gen_on_during_charge==1); %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
    end
    min_on_t = min(min_c,[],2); 
    max_on_t = max(max_c,[],2);
    timedivide_min = max(1,sum(dt(min_on_t~=0)));
    timedivide_max = sum(dt(max_on_t>0));
    min_on_t(isnan(min_on_t))=0;
    max_on_t(isnan(max_on_t))=0;
    mc.Min = sum(min_on_t(min_on_t~=0).*dt(min_on_t~=0))/timedivide_min;%make the cost proportional to amount of time at that cost
    mc.Max = sum(max_on_t(max_on_t>0).*dt(max_on_t>0))/timedivide_max;
end
end%Ends function min_max_cost

function m_c = m_cost(gen,set,scale)
states = gen.states(1:nnz(~cellfun('isempty',gen.states(:,end))),end);
m_c = gen.(states{1}).f(end);
if ~isempty(set) 
    j = 1;
    while j<length(states) && gen.(states{j}).ub(end)<=set
        set = set - gen.(states{j}).ub(end);
        j = j+1;
    end
    m_c = gen.(states{j}).f(end) + set*gen.(states{j}).H(end);
end
if m_c<0
    j = 1;
    while j<length(states) && m_c<0
        j = j+1;
        m_c = gen.(states{j}).f(end);
    end
    if m_c<0
        m_c = 1e-3;
    end
end
m_c = m_c*scale;
end%Ends function m_cost