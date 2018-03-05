function Si = StepDispatchForward(Si,Date,Forecast,Solution)
% Records the last dispatch and steps forward by a single step or a full day
global Plant CurrentState Last24hour DateSim
nS = length(Date)-1;
dt = (24*3600) * (Date(2:end) - Date(1:end-1)); % duration of each time segment [seconds]
if strcmp(Plant.optimoptions.method,'Planning')%assumes forecast is perfect
    DateSim = round(864000*(DateSim+Plant.optimoptions.Horizon/24))/864000;%%count forward by length of the horizon, rounded to nearest second
    prevTime = Last24hour.Timestamp;
    Last24hour = GetCurrentData(Last24hour.Timestamp + Plant.optimoptions.Horizon/24);
    if isfield(Plant,'Building') && ~isempty(Plant.Building)
        nB = length(Plant.Building);
        Tzone = zeros(nS+1,nB);
        Twall = zeros(nS+1,nB);
        for i = 1:1:nB
            [Tzone(:,i),Twall(:,i)] = BuildingSimulate(Plant.Building(i),Forecast.Weather.Tdb,Forecast.Weather.RH,dt,Forecast.Building.InternalGains(:,i),Forecast.Building.ExternalGains(:,i),Solution.Buildings.Cooling(:,i),Solution.Buildings.Heating(:,i),Forecast.Building.AirFlow(:,i),Forecast.Building.Damper(:,i),CurrentState.Buildings(1,i),CurrentState.Buildings(2,i));
        end
        CurrentState.Buildings = [Tzone(end,:);Twall(end,:);DateSim*ones(1,nB)];
        Plant.Design.Buildings(Si+1:Si+nS,:) = Tzone(2:end,:);
    end
    CurrentState.Generators = Solution.Dispatch(end,:);
    CurrentState.Lines = Solution.LineFlows(end,:);
    
    Plant.Design.Timestamp(Si+1:Si+nS) = Forecast.Timestamp;
    Plant.Design.GeneratorState(Si+1:Si+nS,:) = Solution.Dispatch(2:end,:);
    Plant.Design.LineFlows(Si+1:Si+nS,:) = Solution.LineFlows;
    F = fieldnames(Last24hour);
    F = F(~strcmp('Timestamp',F));
    for j = 1:1:length(F)
        if isstruct(Last24hour.(F{j}))
            S = fieldnames(Last24hour.(F{j}));
            for i = 1:1:length(S)
                if ~isempty(Last24hour.(F{j}).(S{i}))
                    Last24hour.(F{j}).(S{i}) = interp1([prevTime;Forecast.Timestamp],[Last24hour.(F{j}).(S{i}); Forecast.(F{j}).(S{i})],Last24hour.Timestamp);
                    Plant.Design.(F{j}).(S{i})(Si+1:Si+nS,:) = Forecast.(F{j}).(S{i});
                end
            end
        elseif ~isempty(Last24hour.(F{j}))
            Last24hour.(F{j}) = interp1([prevTime;Forecast.Timestamp],[Last24hour.(F{j}); Forecast.(F{j})],Last24hour.Timestamp);
            Plant.Design.(F{j})(Si+1:Si+nS,:) = Forecast.(F{j});
        end
    end
    if isfield(Solution,'LBRelax')
        Plant.Design.LBRelax(Si+1:Si+nS) = Solution.LBRelax;
    end
    if isfield(CurrentState,'Hydro')
        CurrentState.Hydro = Solution.hydroSOC(end,:);
        Plant.Design.hydroSOC(Si+1:Si+nS,:) = Solution.hydroSOC;
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Last24hour.Hydro.OutFlow(:,n) = interp1([prevTime;Forecast.Timestamp],[Last24hour.Hydro.OutFlow(:,n); Solution.LineFlows(:,Plant.subNet.Hydro.lineNumber(n))],Last24hour.Timestamp);
            Last24hour.Hydro.SourceSink(:,n) = getHydroFlows(Last24hour.Timestamp,n,'SourceSink');%actual SourceSink according to historical records
        end
    end
    Si = Si + nS; %if 24 hours, take 24 steps
elseif strcmp(Plant.optimoptions.method,'Dispatch') || strcmp(Plant.optimoptions.method,'Control')   
    Si = Si+1;
    DateSim = round(864000*(DateSim+Plant.optimoptions.Resolution/24))/864000;%% count forward 1 step, rounded to nearest second
    Data = GetCurrentData(DateSim);
    Last24hour.Timestamp = [Last24hour.Timestamp(2:end);DateSim;];
    CurrentState.Generators = Solution.Dispatch(2,:);
    CurrentState.Lines = Solution.LineFlows(2,:);
    
    Plant.Dispatch.Timestamp(Si) = DateSim;
    Plant.Dispatch.GeneratorState(Si,:) = CurrentState.Generators;
    Plant.Dispatch.LineFlows(Si,:) = CurrentState.Lines;
    Plant.Predicted.Timestamp(:,Si-1) = Forecast.Timestamp;
    Plant.Predicted.GenDisp(:,:,Si-1) = Solution.Dispatch(2:end,:);
    Plant.Predicted.LineFlows(:,:,Si-1) = Solution.LineFlows;
    Plant.Predicted.Buildings(:,:,Si-1) = Solution.Buildings.Temperature;
    Plant.Predicted.hydroSOC(:,:,Si-1) = Solution.hydroSOC;
    F = fieldnames(Data);
    F = F(~strcmp('Timestamp',F));
    for j = 1:1:length(F)
        if isstruct(Data.(F{j}))
            S = fieldnames(Data.(F{j}));
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
    if isfield(Solution,'LBRelax')
        Plant.Predicted.LBRelax(Si) = Solution.LBRelax;
    end
    if isfield(Plant,'Building') && ~isempty(Plant.Building)
        nB = length(Plant.Building);
        Tzone = zeros(2,nB);
        Twall = zeros(2,nB);
        for i = 1:1:nB
            [Tzone(:,i),Twall(:,i)] = BuildingSimulate(Plant.Building(i),Data.Weather.Tdb(1),Data.Weather.RH(1),dt(1),Data.Building.InternalGains(1,i),Forecast.Building.ExternalGains(1,i),Solution.Buildings.Cooling(1,i),Solution.Buildings.Heating(1,i),Forecast.Building.AirFlow(1,i),Forecast.Building.Damper(1,i),CurrentState.Buildings(1,i),CurrentState.Buildings(2,i));
        end
        CurrentState.Buildings = [Tzone(2,:);Twall(2,:);DateSim*ones(1,nB)];
        Plant.Dispatch.Buildings(Si,:) = CurrentState.Buildings(1,:);
    end
    if isfield(CurrentState,'Hydro')
        CurrentState.Hydro = Solution.hydroSOC(1,:);
        Plant.Dispatch.hydroSOC(Si,:) = CurrentState.Hydro;
        Plant.Predicted.hydroSOC(:,:,Si) = Solution.hydroSOC;
        for n = 1:1:length(Plant.subNet.Hydro.nodes)
            Plant.Dispatch.OutFlow(Si,:) = Solution.LineFlows(1,Plant.subNet.Hydro.lineNumber(n));
            Last24hour.Hydro.OutFlow(:,n) = [Last24hour.Hydro.OutFlow(2:end,n);Solution.LineFlows(1,Plant.subNet.Hydro.lineNumber(n))];
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
end%Ends function StepDispatchForward