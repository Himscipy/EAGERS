function riverFlow = getHydroFlows(Date,node,prop)
%Date is a vector of timstamps that you request the historical inflow rate at.
global Plant
Time = [];
Flow = [];
if all(Date<=Plant.Data.Hydro.Timestamp(1))
    Time = [0 Plant.Data.Hydro.Timestamp(1)];  
    Flow = [Plant.Data.Hydro.(prop)(1,node), Plant.Data.Hydro.(prop)(1,node)]; 
elseif any(Date<=Plant.Data.Hydro.Timestamp(1))
    Time = [0];  
    Flow = [Plant.Data.Hydro.(prop)(1,node)];  
end
if any(Date>Plant.Dispatch.Timestamp(1)) %historical data
    x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
    x2 = nnz(Plant.Data.Hydro.Timestamp<Date(end))+1;
    n = x2-x1+1;
    Time(end+1:end+n) = Plant.Data.Hydro.Timestamp(x1:x2); 
    Flow(end+1:end+n) = Plant.Data.Hydro.(prop)(x1:x2,node);
end
riverFlow = interp1(Time,Flow,Date);
end%Ends function getHydroFlows