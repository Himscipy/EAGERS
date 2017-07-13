function Data = GetHistoricalData(Date)
global Plant 
x1 = max(1,nnz(Plant.Data.Timestamp<Date(1))); 
S = fieldnames(Plant.Data.Demand);
Data.Timestamp = Date;
if length(Date) == 1
    r = (Date - Plant.Data.Timestamp(x1))/(Plant.Data.Timestamp(x1+1) - Plant.Data.Timestamp(x1));
    Data.Temperature = (1-r)*Plant.Data.Temperature(x1) + r*Plant.Data.Temperature(x1+1);
    for i = 1:1:length(S)
        Data.Demand.(S{i}) = (1-r)*Plant.Data.Demand.(S{i})(x1,:) + r*Plant.Data.Demand.(S{i})(x1+1,:);
    end
else
    x2 = max(1,nnz(Plant.Data.Timestamp<Date(end))); 
    Data.Temperature = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Temperature(x1:x2+1),Date);
    for i = 1:1:length(S)
        Data.Demand.(S{i}) = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Demand.(S{i})(x1:x2+1,:),Date);
    end
end
if isfield(Plant.Data,'Hydro')
    %spill, outflow, inflow, and source/sink by node
    if length(Date) == 1
        Data.Hydro.SpillFlow = (1-r)*Plant.Data.Hydro.SpillFlow(x1,:) + r*Plant.Data.Hydro.SpillFlow(x1+1,:);
        Data.Hydro.OutFlow = (1-r)*Plant.Data.Hydro.OutFlow(x1,:) + r*Plant.Data.Hydro.OutFlow(x1+1,:);
        Data.Hydro.InFlow = (1-r)*Plant.Data.Hydro.InFlow(x1,:) + r*Plant.Data.Hydro.InFlow(x1+1,:);
        Data.Hydro.SourceSink = (1-r)*Plant.Data.Hydro.SourceSink(x1,:) + r*Plant.Data.Hydro.SourceSink(x1+1,:);
    else
        Data.Hydro.SpillFlow = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Hydro.SpillFlow(x1:x2+1,:),Date);
        Data.Hydro.OutFlow = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Hydro.OutFlow(x1:x2+1,:),Date);
        Data.Hydro.InFlow = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Hydro.InFlow(x1:x2+1,:),Date);
        Data.Hydro.SourceSink = interp1(Plant.Data.Timestamp(x1:x2+1),Plant.Data.Hydro.SourceSink(x1:x2+1,:),Date);
    end
end
