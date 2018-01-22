function CoSimulation(Date, handles)
global Plant DispatchWaitbar DateSim Last24hour RealTimeData CurrentState
%% Start EnergyPlus cosimulation
installMlep
ep = mlepProcess;
ep.arguments = {'B4', 'USA_VA_Sterling-Washington.Dulles.Intl.AP.724030_TMY3'};
ep.acceptTimeout = 6000;
[status, msg] = ep.start;  
if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end
[status, msg] = ep.acceptSocket;
if status ~= 0
    error('Could not connect to EnergyPlus: %s.', msg);
end

%% EAGERS initialize variables
nG = length(Plant.Generator);%skip initialization
IC = zeros(1,nG);
for i=1:1:nG
    if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        IC(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
end
CurrentState.Generators=IC(1:nG);
CurrentState.Buildings = 22.22;
DateSim = Date;
nS = 365*24/Plant.optimoptions.Resolution+1;
RealTimeData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),nS)';

nS = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
dems = fieldnames(Plant.Data.Demand);
Plant.Dispatch.Temperature = zeros(NumSteps,1);
Plant.Dispatch.Timestamp = zeros(NumSteps,1);
Plant.Dispatch.GeneratorState = zeros(NumSteps,nG);
Plant.Predicted.GenDisp = zeros(nS+1,nG,NumSteps);
Plant.Predicted.Timestamp = zeros(nS+1,NumSteps);
Plant.Predicted.Cost = zeros(NumSteps,1);
Plant.Predicted.Demand = [];
for i = 1:1:length(dems)
    loads = length(RealTimeData.Demand.(dems{i})(1,:));
    Plant.Predicted.Demand.(dems{i}) = zeros(NumSteps,nS,loads);
    Plant.Dispatch.Demand.(dems{i}) = zeros(NumSteps,loads);
    Plant.RunData.Demand.(dems{i}) = zeros(NumSteps,loads);
end

Last24hour = [];
TimeYesterday = linspace(DateSim-1+Plant.optimoptions.Resolution/24,DateSim,24/Plant.optimoptions.Resolution)';
if RealTimeData.Timestamp(1)<(DateSim-1) && RealTimeData.Timestamp(end)>=DateSim
    Last24hour = GetCurrentData(TimeYesterday);
else %need to have this in terms of the first timestep
    Last24hour = GetCurrentData(TimeYesterday+1);
    Last24hour.Timestamp = TimeYesterday;
end
Data = GetCurrentData(DateSim);

Plant.Dispatch.GeneratorState(1,:) = IC;
Plant.Dispatch.Timestamp(1) = DateSim;
Si=1; %counter for # of times dispatch loop has run
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
while Si<NumSteps-1
    %% Read from EnergyPlus
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Parse it to obtain building outputs
    [flag, eptime, outputs] = mlepDecodePacket(packet);
    if flag ~= 0, break; end

    %% run optimization
    Date = DateSim+[0;Time/24];
    if isempty(Plant.Building)
        Data = GetCurrentData(DateSim);
    else
        Data = updateForecast(DateSim,[]);
    end
    Forecast = updateForecast(Date(2:end),Data);%% function that creates demand vector with time intervals coresponding to those selected
    scaleCost = updateGeneratorCost(Date(2:end)); %% All feedstock costs were assumed to be 1 when building matrices 
    [Dispatch,OptSchedule] = NRELoptimization2(CurrentState.Generators,CurrentState.Buildings,Forecast,scaleCost);
    %% record predictions
    for i = 1:1:length(dems)
        loads = length(Forecast.Demand.(dems{i})(1,:));
        for j = 1:1:length(loads)
            Plant.Predicted.Demand.(dems{i})(Si,:,j) = Forecast.Demand.(dems{i})(:,j)';
        end
    end
    SP = [OptSchedule.T_ref(1),OptSchedule.T_ref(1),OptSchedule.Pfc_h(1),Dispatch(2,3)>0];
    %% Write to inputs of E+
    ep.write(mlepEncodeRealData(2, 0, (Si-1)*3600, SP));  
    
    %% Plot to GUI
    if strcmp(get(handles.uipanelMain1,'Visible'),'on')
        backSteps = min(Si,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        History = Plant.Dispatch.GeneratorState(Si-backSteps+1:Si,:);
        HistoryTime = Plant.Dispatch.Timestamp(Si-backSteps+1:Si,:);
        updateGUIstatus(handles,Dispatch(1:2,:),History)
        plotNREL(handles,Date(2:end),Dispatch(2:end,:),HistoryTime,History)
    elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
        Plant.Market.MarginCost = MarginalCapacityCost(Dispatch,Date);
        plotMarginalCapacityCost(handles)
    end
    
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end 
    %% Update status
    CurrentState.Generators = Dispatch(2,1:nG);
    CurrentState.Buildings = 22;
    Plant.Predicted.GenDisp(:,:,Si) = Dispatch;
    Plant.Predicted.Timestamp(:,Si) = Date;
    Si = Si+1;    
    DateSim = round(1e5*(DateSim+Plant.optimoptions.Resolution/24))/1e5;%% count forward 1 step, rounded to nearest second
    Plant.Dispatch.Timestamp(Si) = DateSim;
    Plant.Dispatch.GeneratorState(Si,:) = Dispatch(2,:);
    Last24hour.Timestamp = [Last24hour.Timestamp(2:end);DateSim;];
    F = fieldnames(Last24hour);
    F = F(~strcmp('Timestamp',F));
    for j = 1:1:length(F)
        if isstruct(Last24hour.(F{j}))
            S = fieldnames(Last24hour.(F{j}));
            for i = 1:1:length(S)
                Last24hour.(F{j}).(S{i}) = [Last24hour.(F{j}).(S{i})(2:end,:);Data.(F{j}).(S{i})(1,:);];
                Plant.Dispatch.(F{j}).(S{i})(Si,:) = Forecast.(F{j}).(S{i})(1,:);
            end
        else
            Last24hour.(F{j}) = [Last24hour.(F{j})(2:end,:);Data.(F{j})(1,:);];
            Plant.Dispatch.(F{j})(Si,:) = Forecast.(F{j})(1,:);
        end
    end   
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
end
% Stop EnergyPlus
ep.stop;
disp(['Stopped with flag ' num2str(flag)]);
% ==========FLAGS==============
% Flag	Description
% +1	Simulation reached end time.
% 0	    Normal operation.
% -1	Simulation terminated due to an unspecified error.
% -10	Simulation terminated due to an error during the initialization.
% -20	Simulation terminated due to an error during the time integration.
end%Ends function CoSimulation