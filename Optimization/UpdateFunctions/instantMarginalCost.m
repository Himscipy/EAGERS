function [marginal,StorPower] = instantMarginalCost(FirstProfile,scaleCost,netDemand,IC,dt,t)
global Plant
networkNames = fieldnames(Plant.subNet);
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
All = Out;
nG = length(Plant.Generator);  
marginCost = zeros(1,nG);
CHP = [];
I = zeros(1,nG);
StorPower = zeros(1,nG);
if isfield(netDemand,'H')
    excessHeat = -1*netDemand.H(t);
else
    excessHeat = 0;
end
for i = 1:1:nG
    if Plant.Generator(i).Enabled && ~isempty(Plant.Generator(i).QPform.states)%only use enabled gens
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        if isfield(Plant.Generator(i).QPform,'Stor')
            if ~isempty(FirstProfile)
                loss = dt(t)*(Plant.Generator(i).QPform.Stor.SelfDischarge*Plant.Generator(i).QPform.Stor.UsableSize);
                d_SOC = FirstProfile(t+1,i) - IC(i);%expected output of storage in kWh to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
                if (d_SOC/dt(t) + loss)<0 %discharging
                    StorPower(i) = -(d_SOC/dt(t) + loss)*Plant.Generator(i).QPform.Stor.DischEff;
                else %charging
                    StorPower(i) = -(d_SOC/dt(t) +loss)/Plant.Generator(i).QPform.Stor.ChargeEff; 
                end
            end
            if strcmp(Plant.Generator(i).Source,'Electricity')
                storType.E(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                storType.H(end+1) = i;
                excessHeat = excessHeat + StorPower(i);
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                storType.C(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Water')
                storType.W(end+1) = i;
            end
        elseif ~isempty(states)
            for j = 1:1:length(states)
                if ~isempty(FirstProfile)
                    I(i) = I(i) + min(FirstProfile(t+1,i)-I(i),Plant.Generator(i).QPform.(states{j}).ub(end));
                else I(i) = I(i) + Plant.Generator(i).QPform.(states{j}).ub(end);
                end
            end
            marginCost(i) = Plant.Generator(i).QPform.(states{1}).f(end);
            if ~isempty(FirstProfile) %&& isfield(Plant.Generator(i).QPform,'constCost') %all of these cost terms need to be scaled later on  
                j = 1;
                set = FirstProfile(t+1,i);
                while j<length(states) && Plant.Generator(i).QPform.(states{j}).ub(end)<=set
                    set = set - Plant.Generator(i).QPform.(states{j}).ub(end);
                    j = j+1;
                end
                marginCost(i) = Plant.Generator(i).QPform.(states{j}).f(end) + set*Plant.Generator(i).QPform.(states{j}).H(end);
            end
            if marginCost(i)<0
                j = 1;
                while j<length(states) && marginCost(i)<0
                    j = j+1;
                    marginCost(i) = Plant.Generator(i).QPform.(states{j}).f(end);
                end
                if marginCost(i)<0
                    marginCost(i) = 1e-3;
                end
            end
            S = fieldnames(Plant.Generator(i).QPform.output);
            All.(S{1})(end+1) = i;
        end
        if I(i)>0
            if strcmp(Plant.Generator(i).Type,'CHP Generator')
                S = {'E'};
                CHP(end+1) = i;
                if ~isempty(FirstProfile)
                    j = 1;
                    D = 0;
                    if FirstProfile(t+1,i)>0
                        excessHeat = excessHeat - Plant.Generator(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
                    end
                    while j<=length(states) && FirstProfile(t+1,i)-D>0
                        x_i = min(FirstProfile(t+1,i)-D,Plant.Generator(i).QPform.(states{j}).ub(end));
                        excessHeat = excessHeat + x_i*Plant.Generator(i).QPform.output.H(min(j,length(Plant.Generator(i).QPform.output.H(:,1))),1);
                        D = D + x_i;
                        j = j+1;
                    end
                end
            end
            if isfield(Out,(S{1}))
                Out.(S{1})(end+1) = i;
            end
        end
    end
end
MarginCost = scaleCost(t,:).*marginCost; %marginal cost

if isfield(storType,'E')
    ThisType = MarginCost(Out.E);
    if ~isempty(CHP)
        if excessHeat<=1e-1
            for j = 1:1:length(CHP)
                [~,k] = ismember(CHP(j),Out.E);
                ThisType(k) = 0.6*ThisType(k); %assign 40% of the generator cost to the heat production
            end
        end
    end
    if ~isempty(Out.E)
        marginal.E = max(ThisType);
    else marginal.E = min(MarginCost(All.E));
    end
end

if isfield(storType,'H')
    ThisType = MarginCost(Out.H);
    if ~isempty(CHP)
        ThisType(end+1:end+length(CHP)) = MarginCost(:,CHP)*.4; 
    end
    if ~isempty(ThisType)
        marginal.H = max(ThisType);
    else marginal.H = min(MarginCost(All.H));
    end
end

if isfield(storType,'C')
    Elec = [];
    Absorb = [];
    for j = 1:1:length(All.C)
        Cratio = Plant.Generator(All.C(j)).Output.Cooling(end);
        if strcmp(Plant.Generator(All.C(j)).Source,'Electricity') 
            MarginCost(All.C(j)) = marginal.E/Cratio;
            Elec(end+1) = All.C(j);
        else
            MarginCost(All.C(j)) = marginal.H/Cratio;
            Absorb(end+1) = All.C(j);
        end
    end
    if ~isempty(Out.C)
        marginal.C = max(MarginCost(Out.C));
    else
        marginal.C = min(MarginCost(All.C));
    end
end

if isfield(storType,'W')
    if ~isempty(All.E)
        marginal.W = marginal.E;
    else
        marginal.W = 1;
    end
end

end%Ends function instantMarginalCost