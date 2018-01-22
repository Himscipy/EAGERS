function RunSimulation(Date,handles)
global Model_dir Plant DispatchWaitbar DateSim CurrentState GenAvailTime RestartTime RealTimeData
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
Compare = true;%true  %compare mixed integer to non-mixed integer
nG = length(Plant.Generator);
nB = length(Plant.Building);
DateSim = Date;
loadTestData
WYHydroForecast(Date); %Water Year Forecast for Hydro
if isempty(Plant.Building)
    if any(strcmp(Plant.optimoptions.method,{'Dispatch';'Planning'}))
        timestamp = (DateSim:Plant.optimoptions.Resolution/24:(DateSim+Plant.optimoptions.Interval+Plant.optimoptions.Horizon/24))';
        RealTimeData = interpolateData(timestamp,0.00);%create test data at correct frequency
    else %actual control
        timestamp = (DateSim:Plant.optimoptions.Tmpc/24/3600:(DateSim+Plant.optimoptions.Interval+Plant.optimoptions.Horizon/24))';
        RealTimeData = interpolateData(timestamp,0.00);%create test data at correct frequency
    end
else
    nS = 365*24/Plant.optimoptions.Resolution+1;
    RealTimeData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),nS)';
end

if any(strcmp(Plant.optimoptions.method,{'Control'}))
    for i = 1:1:nG
        if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
            RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
        else
            RestartTime(i) = 0;
        end
    end
    GenAvailTime = ones(1,nG).*DateSim; %  Global variable needed in controller mode
end
if isempty(Plant.OneStep)
    nL = 0;
else
    [~,n] = size(Plant.OneStep.organize);
    nL = n-nG-nB;
end
preAllocateSpace('Dispatch');
nS = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
timers = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
if Compare%%If running for comparisson of mixed and non-mixed integer
    Plant.Predicted.GenDispcQP = zeros(nS,nG+nL+nB,NumSteps);
    Plant.Predicted.CostcQP = zeros(NumSteps,1);
    timers_cQP = zeros(NumSteps,3);
end

if isempty(Plant.Building)
    Data = GetCurrentData(DateSim);
else
    Data = updateForecast(DateSim,[]);
end
resetLast24Hour

% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;

if K ==1
    manualInitialCondition;
else
    automaticInitialCondition(Data);
end
Plant.Dispatch.GeneratorState(1,:) = CurrentState.Generators;
Plant.Dispatch.Timestamp(1) = DateSim;
Si=1; %counter for # of times dispatch loop has run
Solution.Dispatch = [];
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
while Si<NumSteps-1
    
    Date = DateSim+[0;Time/24];
    if Si>1
        D = datevec(Date(end));
        if (D(2) == 10) && (D(3) == 1) && (D(4)<Plant.optmoptions.Resolution) %if october 1st,Run a yearly forecast for hydrology
            WYHydroForecast(Date(end));
        end 
    end 
    
    if isempty(Plant.Building)
        Data = GetCurrentData(DateSim);
    else
        Data = updateForecast(DateSim,[]);
    end
    Forecast = updateForecast(Date(2:end),Data);%% function that creates demand vector with time intervals coresponding to those selected

    if Compare%%If running for comparisson of mixed and non-mixed integer
        Plant.optimoptions.MixedInteger = false;
        Solution_cQP = DispatchLoop(Date,Forecast,Solution.Dispatch);
        timers_cQP(Si,:) = Solution_cQP.timers;
        Plant.optimoptions.MixedInteger = true;
        Plant.Predicted.GenDispcQP(:,:,Si) = Solution_cQP.Dispatch(2:end,:);
        [C,~,~] = NetCostCalc(Solution_cQP.Dispatch,Date,'Dispatch');
        Plant.Predicted.CostcQP(Si) = sum(C);
    end

    Solution = DispatchLoop(Date,Forecast,Solution.Dispatch);
    timers(Si,:) = Solution.timers;
    if ~isempty(handles)
        if strcmp(get(handles.uipanelMain1,'Visible'),'on')
            backSteps = min(Si,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
            History = Plant.Dispatch.GeneratorState(Si-backSteps+1:Si,:);
            HistoryTime = Plant.Dispatch.Timestamp(Si-backSteps+1:Si,:);
            updateGUIstatus(handles,Solution.Dispatch(1:2,:),History)
            plotDispatch(handles,Date(2:end),Solution.Dispatch(2:end,:),HistoryTime,History)
        elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
            Plant.Market.MarginCost = MarginalCapacityCost(Solution.Dispatch,Date);
            plotMarginalCapacityCost(handles)
        end
    end
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end
    [C,~,~] = NetCostCalc(Solution.Dispatch,Date,'Dispatch');
    Plant.Predicted.Cost(Si) = sum(C);
    Si = StepDispatchForward(Si,Date,Forecast,Solution);
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
%     if Compare && rem(Si,100)==0 %%If running for comparisson of mixed and non-mixed integer
%         save(fullfile(Model_dir,'GUI','Optimization','Results',strcat(Plant.Name,'_Compare_MI_',num2str(ceil(Si/100)),'.mat')),'Plant');
%     end
end
% HeatRecovery = CalcHeatRecovery(Plant.Dispatch.GeneratorState(2:end,:));
if strcmp(Plant.optimoptions.method,'Dispatch')
    Plant.Cost.Dispatch = NetCostCalc(Plant.Dispatch.GeneratorState,Plant.Dispatch.Timestamp,'Dispatch');
elseif strcmp(Plant.optimoptions.method,'Control')
    Plant.Cost.RunData = NetCostCalc(Plant.RunData.GeneratorInput,Plant.RunData.Timestamp,'Input');
end