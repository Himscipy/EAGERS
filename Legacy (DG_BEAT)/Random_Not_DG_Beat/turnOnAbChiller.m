function GenOutput = turnOnAbChiller(GenDisp,excessHeat,QP,Locked)
%if there is excess heat available, turn on ab chiller
global Plant
nG = length(Plant.Generator);
nS = length(dt);
abChiller = [];
ColdTank = [];
StoredEnergy = zeros(nS,1);
ChargeEff = [];
StorCapacity = 0;
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(Plant.Generator(i).Source,'Heat')
        abChiller(end+1) = i;
    elseif strcmp(Plant.Generator(i).Type,'Thermal Storage') && strcmp(Plant.Generator(i).Source,'Cooling')
        ColdTank(end+1) = i;
        ChargeEff(end+1) = Plant.Generator(i).QPform.Stor.ChargeEff;
        buff = Plant.Generator(i).QPform.Stor.UsableSize*(Plant.Generator(i).VariableStruct.Buffer/100);
        StoredEnergy = StoredEnergy + GenOutput(2:end,i);
        StorCapacity = StorCapacity + Plant.Generator(i).Size - buff;
    end
end
if ~isempty(abChiller) && ~isempty(ColdTank)
    %%identify times when abChiller could turn on
    avgChargeEff = mean(ChargeEff);
    TurnOn = false(nS,nG);
    NewOut = zeros(nS,nG);
    if isfield(Forecast.Demand,'H')
        excessHeat = -Forecast.Demand.H;
    else
        excessHeat = zeros(nS,1);
    end
    for t = 1:1:nS
        if any(GenOutput(t+1,abChiller)==0) %at least 1 absorption chiller is off
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,{'CHP Generator';'Heater'}) && GenOutput(t+1,i)>0
                    states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
                    D = 0;
                    j = 1;
                    if isfield(Plant.Generator(i).QPform,'constDemand')
                        excessHeat(t) = excessHeat(t) - Plant.Generator(i).QPform.constDemand.H;%for CHP generators, this constDemand is negative
                    end
                    while D<GenOutput(t+1,i) && j<length(states)
                        g = min(GenOutput(t+1,i)-D,Plant.Generator(i).QPform.(states{j}).ub(2));
                        D = D+g;
                        excessHeat(t) = excessHeat(t) + g*Plant.Generator(i).QPform.output.H(j,2);
                        j = j+1;
                    end
                end
            end
            %% estimate possible production from ab chiller at those times
            [~,I] = sort(GenOutput(t+1,abChiller),2,'descend');
            abChiller = abChiller(I);
            for k = 1:1:length(abChiller)
                if GenOutput(t+1,abChiller(k))==0
                    i = abChiller(k);
                    states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
                    for j=1:1:length(states)
                        if j == 1
                            if excessHeat(t)>(Plant.Generator(i).QPform.constDemand.H+Plant.Generator(i).QPform.(states{1}).lb(2)*(-Plant.Generator(i).QPform.output.H(1,2)))%enough heat to get above LB
                                excessHeat(t) = excessHeat(t) - Plant.Generator(i).QPform.constDemand.H;
                                TurnOn(t,i) = true;
                            else
                                excessHeat(t) = 0;
                            end
                        end
                        H = min(excessHeat(t),Plant.Generator(i).QPform.(states{j}).ub(2)*(-Plant.Generator(i).QPform.output.H(j,2)));
                        NewOut(t,i) = NewOut(t,i) + H/(-Plant.Generator(i).QPform.output.H(j,2));
                        excessHeat(t) = excessHeat(t) - H;
                    end
                end
            end
        end
    end
    nT = length(ColdTank);
    for ab = 1:1:length(abChiller)
        i = abChiller(ab);
        for k = 1:1:nnz(TurnOn(:,i))
            %find TurnOn closes to point that is already on
            On = nonzeros((1:nS)'.*(GenOutput(2:end,i)>0));
            t = min(nonzeros((1:nS)'.*TurnOn(:,i)));
            if ~isempty(On)
                d_best = min(abs(t-On));
                for j = 1:1:length(TurnOn(:,i))
                    if TurnOn(j,i) && min(abs(j-On))<d_best
                        t = j;
                        d_best = min(abs(j-On));
                    end
                end
            end
            %if it turns on
            if all((StorCapacity-StoredEnergy(t:end))>NewOut(t,i))
                StoredEnergy(t:end) = StoredEnergy(t:end) + NewOut(t,i)*avgChargeEff;
                GenOutput(t+1,i) = NewOut(t,i);
                addStor = NewOut(t,i);
                for j = 1:1:nT
                    buff = Plant.Generator(ColdTank(j)).QPform.Stor.UsableSize*(Plant.Generator(ColdTank(j)).VariableStruct.Buffer/100);
                    charge = min(addStor/(nT+1-j),(Plant.Generator(ColdTank(j)).Size - buff)/Plant.Generator(ColdTank(j)).QPform.Stor.ChargeEff);
                    GenOutput(t+1,ColdTank(j)) = GenOutput(t+1,ColdTank(j))-charge;
                    addStor = addStor - charge;
                end
            end
            TurnOn(t,i) = false;
        end
    end
end


end%ends function checkAbChiller