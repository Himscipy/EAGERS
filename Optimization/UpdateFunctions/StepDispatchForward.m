function Si = StepDispatchForward(Si,Date,Forecast,Solution)
% Records the last dispatch and steps forward by a single step or a full day
global Plant CurrentState Last24hour DateSim
nG = length(Plant.Generator);
nB = length(Plant.Building);
if isempty(Plant.OpMatA)
    nL = 0;
else
    nL = length(Plant.OpMatA.Organize.IC) - nG - nB;
end
nS = length(Date)-1;

if strcmp(Plant.optimoptions.method,'Planning')
    CurrentState.Generators = Solution.Dispatch(end,1:nG);
    CurrentState.Lines = Solution.Dispatch(end,nG+1:nG+nL);
    CurrentState.Buildings = Solution.Dispatch(end,nG+nL+1:nG+nL+nB);
    Last24hour = GetCurrentData(Last24hour.Timestamp + Plant.optimoptions.Horizon/24);
    Last24hour.Timestamp = Last24hour.Timestamp + Plant.optimoptions.Horizon/24;
    n_last24 = length(Last24hour.Timestamp);
    Plant.Design.GeneratorState(Si:Si+nS,:) = Solution.Dispatch;
    F = fieldnames(Forecast);
    for j = 1:1:length(F)
        if isstruct(Forecast.(F{j}))
            S = fieldnames(Forecast.(F{j}));
            for i = 1:1:length(S)
                Last24hour.(F{j}).(S{i}) = Forecast.(F{j}).(S{i})(end-n_last24+1:end,:);
                Plant.Design.(F{j}).(S{i})(Si+1:Si+nS,:) = Forecast.(F{j}).(S{i});
            end
        else
            Last24hour.(F{j}) = Forecast.(F{j})(end-n_last24+1:end,:);
            Plant.Design.(F{j})(Si+1:Si+nS,:) = Forecast.(F{j});
        end
    end
    Si = Si + nS; %if 24 hours, take 24 steps
    DateSim = round(1e5*(DateSim+Plant.optimoptions.Horizon/24))/1e5;%%count forward by length of the horizon, rounded to nearest second
    if isfield(CurrentState,'Hydro')
        CurrentState.Hydro = Solution.hydroSOC(end,:);
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Last24hour.Hydro.OutFlow(:,n) = Solution.Dispatch(end-n_last24+1:end,nG+Plant.subNet.Hydro.lineNumber(n));
            Last24hour.Hydro.SourceSink(:,n) = getHydroFlows(Last24hour.Timestamp,n,'SourceSink');%actual SourceSink according to historical records
        end
    end
elseif strcmp(Plant.optimoptions.method,'Dispatch') || strcmp(Plant.optimoptions.method,'Control')   
    Si = Si+1;
    DateSim = round(864000*(DateSim+Plant.optimoptions.Resolution/24))/864000;%% count forward 1 step, rounded to nearest second
    Data = GetCurrentData(DateSim);
    CurrentState.Generators = Solution.Dispatch(2,1:nG);
    CurrentState.Lines = Solution.Dispatch(2,nG+1:nG+nL);
    CurrentState.Buildings = Solution.Dispatch(2,nG+nL+1:nG+nL+nB);
    
    Plant.Dispatch.Timestamp(Si) = DateSim;
    Plant.Dispatch.GeneratorState(Si,:) = Solution.Dispatch(2,1:nG);
    Plant.Dispatch.LineFlows = Solution.Dispatch(2,nG+1:nG+nL);
    Plant.Dispatch.Buildings = Solution.Dispatch(2,nG+nL+1:nG+nL+nB);
    Plant.Predicted.GenDisp(:,:,Si-1) = Solution.Dispatch(2:end,:);
    F = fieldnames(Data);
    for j = 1:1:length(F)
        if isstruct(Forecast.(F{j}))
            S = fieldnames(Forecast.(F{j}));
            for i = 1:1:length(S)
                Plant.Dispatch.(F{j}).(S{i})(Si,:) = Data.(F{j}).(S{i})(1,:);
                Last24hour.(F{j}).(S{i}) = [Last24hour.(F{j}).(S{i})(2:end,:);Data.(F{j}).(S{i})(1,:)];
                for k = 1:1:length(Forecast.(F{j}).(S{i})(1,:))
                    Plant.Predicted.(F{j}).(S{i})(Si-1,:,k) = Forecast.(F{j}).(S{i})(:,k)';
                end
            end
        else
            Plant.Dispatch.(F{j})(Si,:) = Data.(F{j})(1,:);
            Last24hour.(F{j}) = [Last24hour.(F{j})(2:end,:);Data.(F{j})(1,:)];
            for k = 1:1:length(Forecast.(F{j})(1,:))
                Plant.Predicted.(F{j})(Si-1,:,k) = Forecast.(F{j})(:,k)';
            end
        end
    end
    if isfield(CurrentState,'Hydro')
        CurrentState.Hydro = Solution.hydroSOC(1,:);
        Plant.Dispatch.SOC(Si,:) = CurrentState.Hydro;
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Last24hour.Hydro.OutFlow(:,n) = [Last24hour.Hydro.OutFlow(2:end,n);Solution.Dispatch(1,nG+Plant.subNet.Hydro.lineNumber(n))];
            Last24hour.Hydro.SourceSink(:,n) = [Last24hour.Hydro.SourceSink(2:end,n);Forecast.Hydro.SourceSink(1,n)];
        end
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