function Si = StepDispatchForward(Si,Date,Data,Forecast,GenDisp)
% Records the last dispatch and steps forward by a single step or a full day
global Plant CurrentState Last24hour DateSim
nG = length(Plant.Generator);
nB = length(Plant.Building);
nL = length(Plant.OpMatA.Organize.IC) - nG - nB;
nS = length(Date)-1;
if isfield(Data,'Demand')
    Outs =  fieldnames(Data.Demand);
else Outs = [];
end

if strcmp(Plant.optimoptions.method,'Planning')
    CurrentState.Generators = GenDisp(end,1:nG);
    CurrentState.Lines = GenDisp(end,nG+1:nG+nL);
    CurrentState.Buildings = GenDisp(end,nG+nL+1:nG+nL+nB);
    if isfield(CurrentState,'Hydro')
        SOC = calculateHydroSOC(GenDisp,Date);
        CurrentState.Hydro = SOC(end,:);
    end
    Last24hour = GetCurrentData(Last24hour.Timestamp + Plant.optimoptions.Horizon/24);
    Last24hour.Timestamp = Last24hour.Timestamp + Plant.optimoptions.Horizon/24;
    Plant.Design.GeneratorState(Si:Si+nS,:) = GenDisp;
    Plant.Design.Timestamp(Si:Si+nS) = Date;
    Plant.Design.Temperature(Si+1:Si+nS) = Forecast.Temperature;
    for i = 1:1:length(Outs)
        Plant.Design.Demand.(Outs{i})(Si+1:Si+nS,:) = Forecast.Demand.(Outs{i});
    end
    Si = Si + nS; %if 24 hours, take 24 steps
    DateSim = round(1e5*(DateSim+Plant.optimoptions.Horizon/24))/1e5;%%count forward by length of the horizon, rounded to nearest second
    if isfield(CurrentState,'Hydro')
        Plant.Data.HydroHistory.Timestamp(end+1:end+nS) = Date(2:end);
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Plant.Data.HydroHistory.OutFlow(end+1:end+nS,n) = GenDisp(2:end,nG+Plant.subNet.Hydro.lineNumber(n));
            Plant.Data.HydroHistory.SourceSink(end+1:end+nS,n) = getHydroFlows(Date(2:end),n,'SourceSink');%actual SourceSink according to historical records
        end
    end
elseif strcmp(Plant.optimoptions.method,'Dispatch') || strcmp(Plant.optimoptions.method,'Control')
    CurrentState.Generators = GenDisp(2,1:nG);
    CurrentState.Lines = GenDisp(2,nG+1:nG+nL);
    CurrentState.Buildings = GenDisp(2,nG+nL+1:nG+nL+nB);
    if isfield(CurrentState,'Hydro')
        global InOut
        SOC = calculateHydroSOC(GenDisp,Date);
        CurrentState.Hydro = SOC(2,:);
        Plant.Dispatch.TransLoss(Si,:) = Plant.Trans(2,:); %%%%
        Plant.Dispatch.InFlow(Si,:) = InOut.InFlow(1,:);
        Plant.Dispatch.OutFlow(Si,:) = InOut.OutFlow(1,:);
        Plant.Dispatch.SS(Si,:) = InOut.SS(1,:);
    end
    Plant.Predicted.GenDisp(:,:,Si) = GenDisp;
    Plant.Predicted.Timestamp(:,Si) = Date;
    Si = Si+1;
    DateSim = round(1e5*(DateSim+Plant.optimoptions.Resolution/24))/1e5;%% count forward 1 step, rounded to nearest second
    Plant.Dispatch.Timestamp(Si) = DateSim;
    Plant.Dispatch.GeneratorState(Si,:) = GenDisp(2,:);
    Plant.Dispatch.Temperature(Si) = Forecast.Temperature(1);
    Last24hour.Timestamp = [Last24hour.Timestamp(2:end);Data.Timestamp;];
    Last24hour.Temperature = [Last24hour.Temperature(2:end,:);Data.Temperature;];
    for i = 1:1:length(Outs)
        Plant.Dispatch.Demand.(Outs{i})(Si,:) = Forecast.Demand.(Outs{i})(1,:);
        Last24hour.Demand.(Outs{i}) = [Last24hour.Demand.(Outs{i})(2:end,:);Data.Demand.(Outs{i})(1,:);];
    end
    for i = 1:1:nB
        Plant.Dispatch.Building(i).Temperature(Si) = GenDisp(2,nG+nL+i);
    end
    if isfield(CurrentState,'Hydro')
        Plant.Data.HydroHistory.Timestamp(end+1) = DateSim;
        nd = length(Plant.Data.HydroHistory.OutFlow);
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Plant.Data.HydroHistory.OutFlow(nd+1,n) = GenDisp(2,nG+Plant.subNet.Hydro.lineNumber(n));
            Plant.Data.HydroHistory.SourceSink(nd+1,n) = getHydroFlows(DateSim,n,'SourceSink');%actual SourceSink according to historical records
        end
        for i = 1:1:nG
            Plant.Dispatch.PowerFlow(Si,i) = GenDisp(2,i)*Plant.Generator(i).QPform.Stor.Power2Flow;
        end
        Plant.Dispatch.SpillFlow(Si,:) = Plant.Dispatch.OutFlow(Si,:) - Plant.Dispatch.PowerFlow(Si,:);
    end
end
if strcmp(Plant.optimoptions.method,'Control')
    %count forward in time in the control loop
        %% Real-time control 
    if any(strcmp(Plant.optimoptions.mode,{'observer';'controller';}))
        if isempty(timerfindall)
            Timers(Plant.optimoptions)
        else %do nothing;
        end
    else
        %% Virtual Plant
        D = DateSim;
        while DateSim<(D+Plant.optimoptions.Resolution/24)
            OnlineLoop
        end
    end
end
% NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution;