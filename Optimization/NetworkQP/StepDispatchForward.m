function Si = StepDispatchForward(Si,Date,Data,Forecast,GenDisp)
% Records the last dispatch and steps forward by a single step or a full day
global Plant CurrentState Last24hour DateSim
nG = length(Plant.Generator);
nS = length(Date)-1;
Outs =  fieldnames(Data.Demand);

if strcmp(Plant.optimoptions.method,'Planning')
    CurrentState.Generators = GenDisp(end,1:nG);
    CurrentState.Lines = GenDisp(end,nG+1:end);
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
elseif strcmp(Plant.optimoptions.method,'Dispatch') || strcmp(Plant.optimoptions.method,'Control')
    CurrentState.Generators = GenDisp(2,1:nG);
    CurrentState.Lines = GenDisp(2,nG+1:end);
    Plant.Predicted.GenDisp(:,:,Si) = GenDisp;
    Plant.Predicted.Timestamp(:,Si) = Date;
    Plant.Dispatch.GeneratorState(Si,:) = [CurrentState.Generators, CurrentState.Lines];
    Plant.Dispatch.Timestamp(Si) = DateSim;
    Plant.Dispatch.Temperature(Si+1) = Forecast.Temperature(1);
    Last24hour.Timestamp = [Last24hour.Timestamp(2:end);Data.Timestamp;];
    Last24hour.Temperature = [Last24hour.Temperature(2:end,:);Data.Temperature;];
    for i = 1:1:length(Outs)
        Plant.Dispatch.Demand.(Outs{i})(Si+1,:) = Forecast.Demand.(Outs{i})(1,:);
        Last24hour.Demand.(Outs{i}) = [Last24hour.Demand.(Outs{i})(2:end,:);Data.Demand.(Outs{i})(1,:);];
    end
    if isfield(Last24hour,'Hydro')
        Last24hour.Hydro.SpillFlow = [Last24hour.Hydro.SpillFlow(2:end,:);Data.Hydro.SpillFlow(1,:);];
        Last24hour.Hydro.OutFlow = [Last24hour.Hydro.OutFlow(2:end,:);Data.Hydro.OutFlow(1,:);];
        Last24hour.Hydro.InFlow = [Last24hour.Hydro.InFlow(2:end,:);Data.Hydro.InFlow(1,:);];
        Last24hour.Hydro.SourceSink = [Last24hour.Hydro.SourceSink(2:end,:);Data.Hydro.SourceSink(1,:);];
    end
    Si = Si+1;
    DateSim = round(1e5*(DateSim+Plant.optimoptions.Resolution/24))/1e5;%% count forward 1 step, rounded to nearest second
    Plant.Dispatch.GeneratorState(Si,:) = GenDisp(2,:);
    Plant.Dispatch.Timestamp(Si) = DateSim;
end
if strcmp(Plant.optimoptions.method,'Control')
    %count forward in time in the control loop
        %% Real-time control 
    if Plant.optimoptions.fastsimulation==0
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