function CoSimulation(date, handles)
global Plant DispatchWaitbar TestData 
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
n_g = length(Plant.Generator);%skip initialization
ic = zeros(1,n_g);
lb = zeros(1,n_g);
for i=1:1:n_g
    switch Plant.Generator(i).Type
        case {'Electric Generator';'CHP Generator';}
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Electricity(end);
        case 'Chiller';
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Cooling(end);
        case 'Heater';
            lb(i) = Plant.Generator(i).VariableStruct.Startup.Heat(end);
        case {'Electric Storage';'Thermal Storage';}
            ic(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
    end
    Plant.Generator(i).CurrentState(1) = ic(i);
    Plant.Generator(i).Status = ic(i)>lb(i);
end
temperatures.build = 22;
n_s = 365*24/Plant.optimoptions.Resolution+1;
TestData.RealTimeData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),n_s)';

n_s = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
num_steps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
dems = fieldnames(Plant.Data.Demand);
Plant.Dispatch.Temperature = zeros(num_steps,1);
Plant.Dispatch.Timestamp = zeros(num_steps,1);
Plant.Dispatch.GeneratorState = zeros(num_steps,n_g);
Plant.Predicted.GenDisp = zeros(n_s+1,n_g,num_steps);
Plant.Predicted.Timestamp = zeros(n_s+1,num_steps);
Plant.Predicted.Cost = zeros(num_steps,1);
Plant.Predicted.Demand = [];
for i = 1:1:length(dems)
    loads = length(TestData.RealTimeData.Demand.(dems{i})(1,:));
    Plant.Predicted.Demand.(dems{i}) = zeros(num_steps,n_s,loads);
    Plant.Dispatch.Demand.(dems{i}) = zeros(num_steps,loads);
    Plant.RunData.Demand.(dems{i}) = zeros(num_steps,loads);
end
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    buildings = Plant.Building;
else
    buildings = [];
end
Plant.Dispatch.GeneratorState(1,:) = ic;
Plant.Dispatch.Timestamp(1) = date;
s_i=1; %counter for # of times dispatch loop has run
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
Time = build_time_vector(Plant.optimoptions);%% set up vector of time interval
while s_i<num_steps-1
    %% Read from EnergyPlus
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Parse it to obtain building outputs
    [flag, eptime, outputs] = mlepDecodePacket(packet);
    if flag ~= 0, break; end

    %% run optimization
    date = date+[0;Time/24];
    
    [forecast,Plant.Generator,buildings] = update_forecast(Plant.Generator,buildings,Plant.subNet,Plant.optimoptions,date(2:end));%% function that creates demand vector with time intervals coresponding to those selected
    scale_cost = update_cost(date(2:end),Plant.Generator); %% All feedstock costs were assumed to be 1 when building matrices 
    [dispatch,opt_schedule] = solver_nrel(Plant.Generator,Plant.optimoptions,GenStatus(1:n_g),temperatures.build,forecast,scale_cost);
    %% record predictions
    for i = 1:1:length(dems)
        loads = length(forecast.Demand.(dems{i})(1,:));
        for j = 1:1:length(loads)
            Plant.Predicted.Demand.(dems{i})(s_i,:,j) = forecast.Demand.(dems{i})(:,j)';
        end
    end
    SP = [opt_schedule.T_ref(1),opt_schedule.T_ref(1),opt_schedule.Pfc_h(1),dispatch(2,3)>0];
    %% Write to inputs of E+
    ep.write(mlepEncodeRealData(2, 0, (s_i-1)*3600, SP));  
    
    %% Plot to GUI
    if strcmp(get(handles.uipanelMain1,'Visible'),'on')
        backSteps = min(s_i,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        history = Plant.Dispatch.GeneratorState(s_i-backSteps+1:s_i,:);
        history_time = Plant.Dispatch.Timestamp(s_i-backSteps+1:s_i,:);
        updateGUIstatus(handles,dispatch(1:2,:),history)
        plotNREL(handles,date(2:end),dispatch(2:end,:),history_time,history)
    elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
        Plant.Market.MarginCost = marginal_cost(Plant.Generator,dispatch,date);
        plotMarginalCapacityCost(handles)
    end
    
    if isempty(DispatchWaitbar)
        return %stop button was pressed
    end 
    %% Update status
    for i=1:1:n_g
        Plant.Generator(i).CurrentState(1) = dispatch(2,1:i);
        Plant.Generator(i).Status = dispatch(2,1:i)>lb(i);
    end
    temperatures.build = opt_schedule.T_ref(1);
    Plant.Predicted.GenDisp(:,:,s_i) = dispatch;
    Plant.Predicted.Timestamp(:,s_i) = date;
    s_i = s_i+1;    
    date = round(1e5*(date+Plant.optimoptions.Resolution/24))/1e5;%% count forward 1 step, rounded to nearest second
    Plant.Dispatch.Timestamp(s_i) = date(1);
    Plant.Dispatch.GeneratorState(s_i,:) = dispatch(2,:);
    F = fieldnames(forecast);
    F = F(~strcmp('Timestamp',F));
    for j = 1:1:length(F)
        if isstruct(forecast.(F{j}))
            S = fieldnames(forecast.(F{j}));
            for i = 1:1:length(S)
                Plant.Dispatch.(F{j}).(S{i})(s_i,:) = forecast.(F{j}).(S{i})(1,:);
            end
        else
            Plant.Dispatch.(F{j})(s_i,:) = forecast.(F{j})(1,:);
        end
    end   
    waitbar(s_i/num_steps,DispatchWaitbar,strcat('Running Dispatch'));
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