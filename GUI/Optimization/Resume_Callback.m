% --- Execute to resume a dispatch
global Model_dir Plant Virtual DispatchWaitbar DateSim mainFig CurrentState 
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)

Virtual = 1;
handles = guihandles(mainFig);
set(handles.Start,'Value',1);%reset start button
set(handles.Stop,'Value',0);%reset stop button
Si = nnz(Plant.Dispatch.Timestamp);
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
dems = fieldnames(Plant.Data.Demand);
weath = {'Tdb';'Twb';'irradDireNorm';};
LastDispatch = [Plant.Predicted.GenDisp(1,:,Si-2);Plant.Predicted.GenDisp(:,:,Si-1);];
CurrentState.Generators = LastDispatch(2,:);
DateSim = Plant.Dispatch.Timestamp(Si);
DispatchWaitbar=waitbar(Si/NumSteps,'Running Dispatch','Visible','off');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
timers = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
if isfield(Plant.Predicted,'GenDispcQP')
    Compare = true; %compare mixed integer to non-mixed integer
    timers_cQP = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
else
    Compare = false;
end


while Si<NumSteps-1
    Date = DateSim+[0;Time/24];
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
% 	if Compare && rem(Si,100)==0 %%If running for comparisson of mixed and non-mixed integer
%         save(fullfile(Model_dir,'GUI','Optimization','Results',strcat(Plant.Name,'_Compare_MI_',num2str(ceil(Si/100)),'.mat')),'Plant');
%     end
end
HeatRecovery = CalcHeatRecovery(Plant.Dispatch.GeneratorState(2:end,:));
if strcmp(Plant.optimoptions.method,'Dispatch')
    Plant.Cost.Dispatch = NetCostCalc(Plant.Dispatch.GeneratorState,Plant.Dispatch.Timestamp,'Dispatch');
elseif strcmp(Plant.optimoptions.method,'Control')
    Plant.Cost.RunData = NetCostCalc(Plant.RunData.GeneratorInput,Plant.RunData.Timestamp,'Input');
end
% Plant.Baseline = RunBaseline(Plant.Dispatch.GeneratorState(1,:)); %finish simulation by running baseline
Stop_Callback(hObject, eventdata, handles)