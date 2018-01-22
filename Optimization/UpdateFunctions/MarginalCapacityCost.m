function MarginCost = MarginalCapacityCost(GenDisp,Date)
%% This function estimates the marginal cost for each addittional kW of generating capacity
%%There are too many factors, e.g. storage, startup, co-production of heat, 
%%reserve capcity called upon at previous step..., to make a complete 
%%determination, so this is purely an estimate. The magin cost is sorted by
%%the cheapest available kW for generators that are already on, then
%%turning on generators one at a time in the cheapest order.
global Plant
n = 4;%break into 4 segments
scaleCost = updateGeneratorCost(Date); 
Binary = GenDisp~=0;
nS = length(Date)-1;
nG = length(Plant.Generator);
nB = length(Plant.Building);
MarginCost.Timestamp = Date;
dt = (Date(2:end) - Date(1:end-1))*24;
S= {'E';'H';'C'};
S2 = {{'CHP Generator';'Electric Generator';};{'Heater'};{'Chiller';}};
S3 = {'Electricity';'Heat';'Cooling';};
MaxOut = GenLimit(GenDisp,Binary,dt);
for k = 1:1:length(S)
    inc = false(nG,1);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,S2{k})
            inc(i) = true;
        end
    end
    if any(inc)
        out = S{k};
        output = S3{k};
        MarginCost.(out).Capacity.SpinReserve = zeros(nG,nS,n);
        MarginCost.(out).Capacity.NonSpin = zeros(nG,nS,n);
        MarginCost.(out).Capacity.DemandResponse = zeros(nG,nS,n);
        MarginCost.(out).Cost.SpinReserve = zeros(nG,nS,n);
        MarginCost.(out).Cost.NonSpin = zeros(nG,nS,n);
        MarginCost.(out).Cost.DemandResponse = zeros(nG,nS,n);
        for i = 1:1:nG
            if inc(i)
                SR = (MaxOut.(S{k})(2:end,i)- GenDisp(2:end,i))';
                for t = 1:1:nS
                    if SR(t)>1e-4
                        Range = linspace(GenDisp(t+1,i),MaxOut.(S{k})(t+1,i),n+1);
                        MarginCost.(out).Capacity.SpinReserve(i,t,:) = SR(t)/n;
                        cost = scaleCost(t+1,i)*Range./interp1(Plant.Generator(i).Output.Capacity*Plant.Generator(i).Size,Plant.Generator(i).Output.(output),Range);
                        MarginCost.(out).Cost.SpinReserve(i,t,:) = (cost(2:n+1)-cost(1:n))/(SR(t)/n);%marginal cost per kW
                        MarginCost.(out).Cost.NonSpin(i,t,:) = inf;
                    else
                        Range = linspace(Plant.Generator(i).VariableStruct.Startup.(output)(end),Plant.Generator(i).Size,n);
                        cost = [0 scaleCost(t+1,i)*Range./interp1(Plant.Generator(i).Output.Capacity*Plant.Generator(i).Size,Plant.Generator(i).Output.(output),Range)];
                        if isfield(Plant.Generator(i).VariableStruct,'StartCost')
                            cost(2) = cost(2) + Plant.Generator(i).VariableStruct.StartCost;
                        end
                        if isfield(Plant.Generator(i).QPform,'constCost')
                            cost(2) = cost(2) + Plant.Generator(i).QPform.constCost*scaleCost(t+1,i);
                        end
                        if strcmp(Plant.Generator(i).Type,'Chiller') && isfield(Plant.Generator(i).QPform.constDemand,'E')
                            cost(2) = cost(2) + Plant.Generator(i).QPform.constDemand.E*min(nonzeros(MarginCost.E.Cost.SpinReserve(:,t,1)));
                        end
                        MarginCost.(out).Capacity.NonSpin(i,t,:) = Range - [0, Range(1:n-1)];
                        MarginCost.(out).Cost.NonSpin(i,t,:) = (cost(2:n+1)-cost(1:n))'./squeeze(MarginCost.(out).Capacity.NonSpin(i,t,:));%marginal cost per kW
                        MarginCost.(out).Cost.SpinReserve(i,t,:) = inf;
                    end
                end
            end
        end

    end 
end
%%make corrections for storage?
    
    %%find cost per kW of demand response
    
end %End function MarginalCapacityCost