function RunPlanning
global testSystems Plant TestData DispatchWaitbar Virtual RealTime DateSim CurrentState
Virtual = 1;
RealTime = 0;
handles = guihandles;
list = get(handles.popupmenuAxes,'String');
item = list{get(handles.popupmenuAxes,'Value')};
if strcmp(item,'Monthly Costs')
    handlesAll = guihandles;
    handlesAll.axesMain = handles.axesMain;
    handlesAll.axesCumulative = handles.axesCumulative;
    handles = handlesAll; %workaround to pass both axes handles and showSys handles
    set(handles.sliderZoom1,'Visible','off');set(handles.sliderDate1,'Visible','off');set(handles.textDate1,'Visible','off');
    set(handles.textDay1,'Visible','off'); set(handles.textAllData1,'Visible','off'); set(handles.textHorizon1,'Visible','off');
    Years = str2double(get(handles.NPC_Years,'String'));%choose years in GUI
    for i_ts = 1:1:length(testSystems)%Run through list of projects
        if ~isfield(testSystems(i_ts),'Design') || ~isfield(testSystems(i_ts).Design,'Timestamp') || any(testSystems(i_ts).Design.Timestamp==0)%if the project has already been run, don't re-simulate (need to empty Design when something in the GUI changes what the solution should be
            Plant = testSystems(i_ts);
            Plant.optimoptions.method = 'Planning';
            Plant.optimoptions.forecast = 'Perfect';%Perfect forecast pulls directly from TestData
            if isfield(Plant,'Building') && ~isempty(Plant.Building)
                Plant.optimoptions.forecast = 'Building';
            end
            Plant.optimoptions.Interval = floor(TestData.Timestamp(end)-TestData.Timestamp(1));
            if get(handles.DesignDay,'Value') ==1%If design days option is selected, optimize the 1st day of the month, and assume the rest of the month to be identical
                Plant.optimoptions.endSOC = 'Initial';%Constrain the final SOC of any storage device to be equal to the initial charge so that days 2-30 of each month do not over depleate storage
                Plant.optimoptions.Horizon = max(24,Plant.optimoptions.Horizon);%Make the horizon at least 1 day
                Plant.subNet = [];%empty the optimization matrices so that they are rebuilt with the new endSOC constraint
                DateSim = TestData.Timestamp(1);%set the starting date
                initializeOptimization%load optimizations
                NumSteps = preAllocateSpace('Design');%create Plant.Design structure with correct space
                interpolateData(Plant.optimoptions.Resolution*3600,Plant.optimoptions.Interval,0.00);%create test data at correct frequency
                STR = 'Optimizing Design Day Dispatch';
                DispatchWaitbar=waitbar(0,STR,'Visible','on');
                reloadLast24hour(DateSim,Plant.optimoptions.Resolution)%re-load the previous 24 hours
                automaticInitialCondition(GetCurrentData(DateSim)); 
                Plant.Design.Timestamp(1) = DateSim;
                ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
                Si = 1;
                while DateSim+Plant.optimoptions.Horizon/24<=TestData.Timestamp(end)%loop to simulate days 1 to n in TestData
                    D = datevec((DateSim));
                    if Si == 1 || D(3) == 1  %If it is the first step, or the first of the month run the actual optimization                     
                        Forecast = updateForecast(ForecastTime(2:end));%% function that creates demand vector with time intervals coresponding to those selected
                        Solution = DispatchLoop(ForecastTime,Forecast,[]);
                    else%otherwise just change the dates and use the previous solution,
                        Forecast.Timestamp = ForecastTime(2:end);
                    end                   
                    StepDispatchForward(Si,ForecastTime,Forecast,Solution);%put solution into Plant.Design
                    if isfield(CurrentState,'Buildings')%don't acumulate error in BuildingSimulate, force warmup
                        CurrentState.Buildings(3,:) = 0;
                    end
                    ForecastTime = round(864000*(ForecastTime+Plant.optimoptions.Horizon/24))/864000;%%count forward by length of the horizon, rounded to nearest second
                    Si = Si + length(ForecastTime)-1;
                    waitbar(Si/NumSteps,DispatchWaitbar,strcat('Running Design Day Dispatch'));
                end
                close(DispatchWaitbar)
                DispatchWaitbar = [];
            else
                Plant.optimoptions.endSOC = 'Flexible';%remove constraint on final SOC of storage
                Plant.subNet = [];%empty the optimization matrices so that they are rebuilt without endSOC constraint
                initializeOptimization
                preAllocateSpace('Design')
                if ~isfield(Plant,'Design') || isempty(Plant.Design) || any(Plant.Design.Timestamp==0) %at least some points have not been run
                    STR = 'Optimizing Dispatch Throughout Entire Year';
                    DispatchWaitbar=waitbar(0,STR,'Visible','on');
                    RunSimulation(TestData.Timestamp(1),[]);
                    close(DispatchWaitbar)
                    DispatchWaitbar = [];
                end
            end
            Plant.optimoptions.method = testSystems(i_ts).optimoptions.method;%change back method so it is correct when switching to control tool
            Plant.optimoptions.forecast = testSystems(i_ts).optimoptions.forecast;%change back method so it is correct when switching to control tool
            testSystems(i_ts).Design = Plant.Design;
        end
        DesignCosts(i_ts,Years,testSystems(i_ts).Costs.Equipment);%update the costs, monthly costs & NPC for system k
    end
    base = get(handles.popupmenuBaseline,'Value');%Which project is the baseline, determined  by GUI
    for i_ts = 1:1:length(testSystems)%calculate caparative economic terms
        if i_ts == base
            testSystems(i_ts).Costs.Financial.irr = 0;
            testSystems(i_ts).Costs.Financial.payback = 0;
        else
            [testSystems(i_ts).Costs.Financial.irr,testSystems(i_ts).Costs.Financial.payback] = FinancialMetric(testSystems(i_ts).Costs.ProjectedMonthlyCosts,testSystems(base).Costs.ProjectedMonthlyCosts,(1+testSystems(i_ts).Costs.DiscountRate/100));
        end
    end
    %plot results
    PlotCosts(handles)
else
    set(handles.sliderZoom1,'Visible','on');set(handles.sliderDate1,'Visible','on');set(handles.textDate1,'Visible','on');
    set(handles.textDay1,'Visible','on'); set(handles.textAllData1,'Visible','on'); set(handles.textHorizon1,'Visible','on');
    PlotHistorical(handles,item,1)
end