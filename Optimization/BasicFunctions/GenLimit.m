function [MaxOut,Constraint,spareGenCumulative] = GenLimit(GenOutput,Binary,dt)
%% Find Max output of each generator at each time
global Plant
[n,nG] = size(Binary);
nS = n - 1;
RampUp = zeros(nS,nG);
RampDown = zeros(nS,nG);
UB = zeros(1,nG);
LB = zeros(1,nG);
Constraint = zeros(nS,nG);
for i = 1:1:nG
    if isfield(Plant.Generator(i).QPform,'Ramp')
        RampUp(:,i) = Plant.Generator(i).QPform.Ramp.b(1)*dt;
        RampDown(:,i) = Plant.Generator(i).QPform.Ramp.b(2)*dt;
        starts = nonzeros((1:nS)'.*((Binary(2:end,i)-Binary(1:nS,i))>0));
        stops = nonzeros((1:nS)'.*(Binary(1:nS,i)-(Binary(2:end,i))>0));
        LB(i) = Plant.Generator(i).QPform.(Plant.Generator(i).QPform.states{1,end}).lb(end);
        RampUp(starts,i) = max(RampUp(starts,i),LB(i)); %if the ramp rate is less than the lb, increase ramp rate at moment of startup
        RampDown(stops,i) = max(RampDown(stops,i),LB(i)); %if the ramp rate is less than the lb, increase ramp rate at moment of shutdown
        UB(i) = Plant.Generator(i).Size;
    end
end
S= {'E';'H';'C'};
S2 = {{'CHP Generator';'Electric Generator';};{'Heater'};{'Chiller';}};
for k = 1:1:length(S)
    spareGen = zeros(nS,1);
    spareGenCumulative.(S{k}) = zeros(nS,1);
    MaxOut.(S{k}) = zeros(nS+1,nG);
    inc = false(nG,1);
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,S2{k})
            inc(i) = true;
        end
    end
    for i = 1:1:nG
        if inc(i)
            MaxOut.(S{k})(1,i) = GenOutput(1,i);
            start = inf; %constrained by initial condition (can't do anything)
            for t = 1:1:nS
                if Binary(t+1,i)
                    if~Binary(t,i)
                        start = t;%just turned on
                    end
                    if UB(i)<=(MaxOut.(S{k})(t,i) + RampUp(t,i))
                        MaxOut.(S{k})(t+1,i) = UB(i);
                    else
                        MaxOut.(S{k})(t+1,i) = (MaxOut.(S{k})(t,i) + RampUp(t,i));
                        if ~isinf(start) && start>1
                            Constraint(t,i) = start-1;
                        else
                            Constraint(t,i) = inf;
                        end
                    end
                end
            end
            for t = nS:-1:2
                if MaxOut.(S{k})(t,i) - RampDown(t,i) > MaxOut.(S{k})(t+1,i) %constrained by shutdown
                    MaxOut.(S{k})(t,i) = min(UB(i),MaxOut.(S{k})(t+1,i) + RampDown(t,i));
                    if Constraint(t-1,i) == 0
                        Constraint(t-1,i) = nS +1 - max(((nS-t+1):-1:1)'.*(MaxOut.(S{k})(t+1:nS+1,i)==0));
                    end
                end
            end
            if strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(Plant.Generator(i).Source,'Heat')%assume net heat production remains the same, so do not count spare absorption chiller capacity
                MaxOut.(S{k})(2:end,i) = GenOutput(2:end,i);%don't count towards spare capacity
            end
            spareGen = spareGen + (MaxOut.(S{k})(2:end,i) - GenOutput(2:end,i));
        end
        if strcmp(S{k},'H') && strcmp(Plant.Generator(i).Type,'CHP Generator')
            heatOutput = zeros(nS+1,1);
            %convert MaxOut.E to MaxOut.H
            %Convert GenOutput to heatOutput
            for t = 0:1:nS
                D = 0;
                H = 0;
                j = 1;
                if MaxOut.E(t+1,i)>0
                    heatOutput(t+1) = -Plant.Generator(i).QPform.constDemand.H;
                    MaxOut.H(t+1,i) = -Plant.Generator(i).QPform.constDemand.H;                
                else
                    heatOutput(t+1) = 0;
                    MaxOut.H(t+1,i) = 0;
                end
                states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
                while j<=length(states)
                    d = min(GenOutput(t+1,i) - D,Plant.Generator(i).QPform.(states{j}).ub(2));
                    D = D + d;
                    heatOutput(t+1) = heatOutput(t+1) + d*Plant.Generator(i).QPform.output.H(j,2);
                    h = min(MaxOut.E(t+1,i) - H,Plant.Generator(i).QPform.(states{j}).ub(2));
                    H = H + h;
                    MaxOut.H(t+1,i) = MaxOut.H(t+1,i) + h*Plant.Generator(i).QPform.output.H(j,2);
                    j = j+1;
                end
            end
            spareGen = spareGen + (MaxOut.(S{k})(2:end,i) - heatOutput(2:end));
        end
    end
    for t = 1:1:nS
        spareGenCumulative.(S{k})(t) = sum(spareGen(1:t).*dt(1:t));
    end
end
end%ends function GenLimit