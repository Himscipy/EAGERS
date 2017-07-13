function Data = GetCurrentData(Date)
global RealTimeData RealTime Virtual
if ~isinf(Date(1)) && (isempty(RealTimeData) || all(Date<=RealTimeData.Timestamp(end)))
    %% need better indexing system to call current data
    x1 = max(1,nnz(RealTimeData.Timestamp<Date(1))); 
    S = fieldnames(RealTimeData.Demand);
    Data.Timestamp = Date;
    if length(Date) == 1
        r = (Date - RealTimeData.Timestamp(x1))/(RealTimeData.Timestamp(x1+1) - RealTimeData.Timestamp(x1));
        Data.Temperature = (1-r)*RealTimeData.Temperature(x1) + r*RealTimeData.Temperature(x1+1);
        for i = 1:1:length(S)
            Data.Demand.(S{i}) = (1-r)*RealTimeData.Demand.(S{i})(x1,:) + r*RealTimeData.Demand.(S{i})(x1+1,:);
        end
    else
        x2 = max(1,nnz(RealTimeData.Timestamp<Date(end))); 
        Data.Temperature = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Temperature(x1:x2+1),Date);
        for i = 1:1:length(S)
            Data.Demand.(S{i}) = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Demand.(S{i})(x1:x2+1,:),Date);
        end
    end
    if isfield(RealTimeData,'Hydro')
        %spill, outflow, inflow, and source/sink by node
        if length(Date) == 1
            Data.Hydro.SpillFlow = (1-r)*RealTimeData.Hydro.SpillFlow(x1,:) + r*RealTimeData.Hydro.SpillFlow(x1+1,:);
            Data.Hydro.OutFlow = (1-r)*RealTimeData.Hydro.OutFlow(x1,:) + r*RealTimeData.Hydro.OutFlow(x1+1,:);
            Data.Hydro.InFlow = (1-r)*RealTimeData.Hydro.InFlow(x1,:) + r*RealTimeData.Hydro.InFlow(x1+1,:);
            Data.Hydro.SourceSink = (1-r)*RealTimeData.Hydro.SourceSink(x1,:) + r*RealTimeData.Hydro.SourceSink(x1+1,:);
        else
            Data.Hydro.SpillFlow = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Hydro.SpillFlow(x1:x2+1,:),Date);
            Data.Hydro.OutFlow = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Hydro.OutFlow(x1:x2+1,:),Date);
            Data.Hydro.InFlow = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Hydro.InFlow(x1:x2+1,:),Date);
            Data.Hydro.SourceSink = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.Hydro.SourceSink(x1:x2+1,:),Date);
        end
    end
else % End Operation or Simulation
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
    Data = [];
end