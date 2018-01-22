function ForecastPlot
global Plant DateSim GENINDEX
handles = guihandles;
%find the current date
if isempty(Plant.Building)
    DateSim = Plant.Data.Timestamp(1) + get(handles.sliderDate,'Value');
    resetLast24Hour
    Data = GetHistoricalData(DateSim);  
else
    DateSim = datenum([2017,1,1]) + get(handles.sliderDate,'Value');
    resetLast24Hour
    Data = updateForecast(DateSim,[]);
end

      
if isempty(GENINDEX)%plot demands
    set(handles.sliderZoom,'Visible','on')
    set(handles.textHour,'Visible','on'); set(handles.textWeek,'Visible','on'); set(handles.textHorizon,'Visible','on');
    nPlot = 0;
    while isfield(handles,strcat('ForecastPlot',num2str(nPlot+1)))
        nPlot = nPlot+1;
        set(handles.(strcat('ForecastPlot',num2str(nPlot))),'Visible','on');%plotting axes (y-axis on left)
        set(handles.(strcat('ForecastPlot',num2str(nPlot),'b')),'Visible','on');%y-axis on the right
        set(handles.(strcat('ForecastName',num2str(nPlot))),'Visible','on');%button with name
    end
    Z = get(handles.sliderZoom,'Value');
    Z_max = get(handles.sliderZoom,'Max');
    if isfield(Plant.Data,'Timestamp')
        lastDate = Plant.Data.Timestamp(end);
    else
        lastDate = datenum([2018,1,1]);
    end
    if Z == Z_max %show all data
        DateEnd = lastDate;
    elseif Z<1
        DateEnd = min(lastDate,DateSim + 1);
    elseif Z<2
        DateEnd = min(lastDate,DateSim + 7);
    elseif Z<3
        DateEnd = min(lastDate,DateSim + 31);
    elseif Z<4
        DateEnd = min(lastDate,DateSim + 365);
    end
        
    Forecast = updateForecast(linspace(DateSim,DateEnd)',Data);%% function that creates demand vector with time intervals coresponding to those selected
    if isempty(Plant.Building)
        Outs =  fieldnames(Forecast.Demand);
        Xi = nnz(Plant.Data.Timestamp<(DateSim-1))+1;
        Xf = nnz(Plant.Data.Timestamp<(DateEnd) & Plant.Data.Timestamp>0);
        for i = 1:1:length(Outs)
            Actual.(Outs{i}) = Plant.Data.Demand.(Outs{i})(Xi:Xf);
        end
        ActualTime = Plant.Data.Timestamp(Xi:Xf);
    else
        %% need to add heating and cooling
        nB = length(Plant.Building);
        if ~isfield(Forecast,'Demand') || ~isfield(Forecast.Demand,'E')
            Forecast.Demand.E = 0;
        end
        Forecast.Demand.T = zeros(length(Forecast.Building(1).E0),nB);
        for i = 1:1:nB
            Forecast.Demand.E = Forecast.Demand.E + Forecast.Building(i).E0;
            Forecast.Demand.T(:,i) = Forecast.Building(i).Tset;
        end
        Actual = [];
        ActualTime = [];
    end
    plotDemand(Forecast.Timestamp,Forecast.Demand,ActualTime,Actual)
else
    set(handles.sliderZoom,'Value',Plant.optimoptions.Horizon/24/7,'Visible','off')%hide slider
    set(handles.textHour,'Visible','off'); set(handles.textWeek,'Visible','off'); set(handles.textHorizon,'Visible','off');
    if ~isfield(Plant,'Dispatch') || isempty(Plant.Dispatch) || ~isfield(Plant.Dispatch,'Timestamp') || isempty(Plant.Dispatch.Timestamp) || min(abs(Plant.Dispatch.Timestamp-DateSim))>=Plant.optimoptions.Resolution/24
        ForecastTime = DateSim+[0;buildTimeVector(Plant.optimoptions)/24];%linspace(DateSim,DateEnd)';would need to re-do optimization matrices for this time vector
        Solution = SingleOptimization(ForecastTime,[]);
        Dispatch = Solution.Dispatch;
        History = [];
        HistoryTime = [];
    else
        t = max(linspace(1,length(Plant.Dispatch.Timestamp))'.*(Plant.Dispatch.Timestamp<DateSim)); %index preceeding current step
        if Plant.Dispatch.Timestamp(t+1)>0 && (Plant.Dispatch.Timestamp(t+1)-DateSim)<(DateSim-Plant.Dispatch.Timestamp(t))
            t = t+1; %The next time step is actually closer
        end
        ForecastTime = Plant.Predicted.Timestamp(t,:)';
        Dispatch = Plant.Predicted.GenDisp(:,:,t);
        backSteps = min(t,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        History = Plant.Dispatch.GeneratorState(t-backSteps+1:t,:);
        HistoryTime = Plant.Dispatch.Timestamp(t-backSteps+1:t);
    end
    Forecast = Dispatch(:,GENINDEX);
    History = History(:,GENINDEX);
    if isfield(Plant.Data,'Dispatch') && any((Plant.Data.Timestamp>0)&(Plant.Data.Timestamp<=DateSim)) && any((Plant.Data.Timestamp>0)&(Plant.Data.Timestamp>=(DateSim + get(handles.sliderZoom,'Value')*7)))
        Xi = nnz(Plant.Data.Timestamp<DateSim);
        Xf = nnz(Plant.Data.Timestamp<(DateSim + get(handles.sliderZoom,'Value')*7));
        Actual = Plant.Data.Dispatch(Xi:Xf,GENINDEX);
        ActualTime = Plant.Data.Timestamp(Xi:Xf);
    else
        Actual = []; ActualTime = [];
    end
    plotSchedule(ForecastTime,Forecast,HistoryTime,History,ActualTime,Actual,GENINDEX)
end
end %Ends ForecastPlot