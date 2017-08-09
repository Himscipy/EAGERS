function SourceSink = getHydroSourceSink(Date,node)
global Plant
x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
SourceSink = zeros(length(Date),1);
if length(Date) == 1
    r = (Date - Plant.Data.Hydro.Timestamp(x1))/(Plant.Data.Hydro.Timestamp(x1+1) - Plant.Data.Hydro.Timestamp(x1));
    SourceSink = (1-r)*Plant.Data.Hydro.SourcesandSinks(x1,node) + r*Plant.Data.Hydro.SourcesandSinks(x1+1,node);
else
    x2 = nnz(Plant.Data.Hydro.Timestamp<Date(end))+1; 
    n = nnz(Date<=Plant.Data.Hydro.Timestamp(x1));%take care of any initial conditions before data exists
    SourceSink(1:n) = Plant.Data.Hydro.Outflow(x1);
    SourceSink(n+1:end) = interp1(Plant.Data.Hydro.Timestamp(x1:x2+1),Plant.Data.Hydro.SourcesandSinks(x1:x2+1,node),Date(n+1:end));
end