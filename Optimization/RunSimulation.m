function RunSimulation(Date,handles)
global Plant DispatchWaitbar DateSim RealTime Virtual
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
DateSim = Date;
if ~isfield(Plant,'subNet') || isempty(Plant.subNet)
    initializeOptimization
end
if isfield(Plant,'Dispatch') && any(Plant.Dispatch.Timestamp)
    Si = nnz(Plant.Dispatch.Timestamp);
    if Si>1 && Plant.Dispatch.Timestamp(Si) == DateSim %resuming after sucessful hitting stop button
        NumSteps = length(Plant.Dispatch.Timestamp); %number of simulation steps
    elseif Si>1 && Plant.Dispatch.Timestamp(Si-1) == DateSim %resuming after crash, hadn't updated plot yet
        NumSteps = length(Plant.Dispatch.Timestamp); %number of simulation steps
    else
        Si = 1;
        NumSteps = preAllocateSpace('Dispatch');
        startNewSimulation
    end
else
    Si = 1;
    NumSteps = preAllocateSpace('Dispatch');
    startNewSimulation
end
timers = zeros(NumSteps,3); % To record times set to zeros(1,3), to not record set to [];
Solution.Dispatch = [];
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
Time = buildTimeVector(Plant.optimoptions);%% set up vector of time interval
while Si<NumSteps-1
    Date = DateSim+[0;Time/24];
    Forecast = updateForecast(Date(2:end));%% function that creates demand vector with time intervals coresponding to those selected
    Solution = DispatchLoop(Date,Forecast,Solution);
    timers(Si,:) = Solution.timers;
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end
    [C,~,~] = NetCostCalc(Solution.Dispatch,Date,'Dispatch');
    Plant.Predicted.Cost(Si) = sum(C);
    Si = StepDispatchForward(Si,Date,Forecast,Solution);
    updateGUIstatus(handles,Solution,Date,Si-1)
    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Dispatch'));
end
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