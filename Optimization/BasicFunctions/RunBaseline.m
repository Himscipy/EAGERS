function Baseline = RunBaseline(IC)
% the goal here is to run simple dynamic economic dispatch
% we can explore other alternatives for what the baseline is
global Plant
Outs = fieldnames(Plant.Data.Demand);
Date = Plant.Data.Timestamp(1);
options = Plant.optimoptions;
options.tspacing = 'constant';
Time = buildTimeVector(options);
dt = Time - [0; Time(1:end-1)];
nG = length(Plant.Generator);
Demand.Renewable = zeros(length(Time),nG);
t = 1;
steps = nnz(Plant.Data.Timestamp<Plant.Data.Timestamp(1)+Plant.optimoptions.Horizon/24); %number of data points in horizon
Baseline.NetCost = 0;
BaselineWaitbarHandle=waitbar(0,'Running Baseline Dynamic Economic Dispatch');
while Date<Plant.Data.Timestamp(1)+Plant.optimoptions.Interval
    Baseline.Timestamp = Date+[0,Time]./24;
    scaleCost = updateGeneratorCost(Baseline.Timestamp);
    margincost = updateMarginalCost(IC,scaleCost,dt);

    for S = 1:1:length(Outs)
        Demand.(Outs{S}) = Plant.Data.Demand.(Outs{S})(t:t+nS-1);
    end
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Demand.Renewable(:,i) = RenewableOutput(i,Date,Time,'Actual');
        end
    end
    QP = updateMatrices(Plant.OpMatB.QP,IC,Baseline.Timestamp,scaleCost,margincost,Demand,[]); 
    
    Locked = ones(length(Time),nG); % will keep all generators on at all times. May not be feasible
    [GenDisp,~,Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B, and all generators on
    
    if Feasible ==1
        Baseline.GeneratorState(t:t+nS,:) = GenDisp;
        Baseline.NetCost = Baseline.NetCost + NetCostCalc(GenDisp,Baseline.Timestamp,'Dispatch');
        t = t+ steps;
        Date = Date + Plant.optimoptions.Horizon/24;
        IC = GenDisp(end,:);
    else
        disp('baseline dispatch is infeasible with all generators on')
        Date = Plant.Data.Timestamp(1)+Plant.optimoptions.Interval+1;
    end
end
%%--%%
close(BaselineWaitbarHandle)