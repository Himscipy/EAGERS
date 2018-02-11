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
Solution.Dispatch = [Plant.Predicted.GenDisp(1,:,Si-2);Plant.Predicted.GenDisp(:,:,Si-1);];
CurrentState.Generators = Solution.Dispatch(2,:);
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
    if Si>1
        D = datevec(Date(end));
        if (D(2) == 10) && (D(3) == 1) && (D(4)<Plant.optimoptions.Resolution) %if october 1st,Run a yearly forecast for hydrology
            WYHydroForecast(Date(end));
        end 
    end 
    Forecast = updateForecast(Date(2:end));%% function that creates demand vector with time intervals coresponding to those selected
    if Compare%%If running for comparisson of mixed and non-mixed integer
        Plant.optimoptions.MixedInteger = false;
        Solution_cQP = DispatchLoop(Date,Forecast,Solution);
        timers_cQP(Si,:) = Solution_cQP.timers;
        Plant.optimoptions.MixedInteger = true;
        if ~isempty(Plant.Generator)
            Plant.Predicted.GenDispcQP(:,:,Si) = Solution_cQP.Dispatch(2:end,:);
        end
        [C,~,~] = NetCostCalc(Solution_cQP.Dispatch,Date,'Dispatch');
        Plant.Predicted.CostcQP(Si) = sum(C);
    end

    Solution = DispatchLoop(Date,Forecast,Solution);
    timers(Si,:) = Solution.timers;   
    updateGUIstatus(handles,Solution,Date,Si)
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end
    [C,~,~] = NetCostCalc(Solution.Dispatch,Date,'Dispatch');
    Plant.Predicted.Cost(Si) = sum(C);
    Si = StepDispatchForward(Si,Date,Forecast,Solution);
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
end
% HeatRecovery = CalcHeatRecovery(Plant.Dispatch.GeneratorState(2:end,:));
if strcmp(Plant.optimoptions.method,'Dispatch')
    Plant.Cost.Dispatch = NetCostCalc(Plant.Dispatch.GeneratorState,Plant.Dispatch.Timestamp,'Dispatch');
elseif strcmp(Plant.optimoptions.method,'Control')
    Plant.Cost.RunData = NetCostCalc(Plant.RunData.GeneratorInput,Plant.RunData.Timestamp,'Input');
end
if RealTime
    closePorts;
end
RealTime=0;%end condition for real simulation
Virtual = 0;%end condition for virtual simulation

T1 = timerfind('Name', 'dispTimer') ;
T2 = timerfind('Name', 'optTimer') ;
T3 = timerfind('Name', 'mpcTimer') ;
T4 = timerfind('Name', 'fanTimer') ;
Timers = [T1,T2,T3,T4];
for i = 1:1:length(Timers)
    stop(Timers(i));
    delete(Timers(i))
end