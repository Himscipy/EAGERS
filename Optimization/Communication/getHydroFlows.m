function riverFlow = getHydroFlows(Date,line)
%Date is a vector of timstamps that you request the flow rate at. If we are
%at hour 5, and you need the previous 8 hours of flow data, 4 points will
%be from previous dispatchs, and 4 from the historical record. The last
%value in Date should be 1 timestep ago.
% line is the index of hydro line names
global Plant
nG = length(Plant.Generator);
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
name = Plant.subNet.lineNames.Hydro{line};
r = strfind(name,'_');
I = find(strcmp(name(1:r(1)-1),Plant.Data.Hydro.Nodes));
Time = [];
Flow = [];
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        nLcum = nLcum+nLinet(net); 
    else
        if isempty(Plant.Dispatch)
            Time = [0 Date(end)];  
            if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                Flow = [Plant.Data.Hydro.OutFlow(1,I), Plant.Data.Hydro.OutFlow(1,I)];  
            else Flow = [Plant.Data.Hydro.SpillFlow(1,I), Plant.Data.Hydro.SpillFlow(1,I)]; 
            end
        else
            if any(Date<=Plant.Dispatch.Timestamp(1)) %historical data
                if any(Date<=Plant.Data.Hydro.Timestamp(1))
                    Time = 0;
                    if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                        Flow = Plant.Data.Hydro.OutFlow(1,I);  
                    else Flow = Plant.Data.Hydro.SpillFlow(1,I); 
                    end
                end
                
                x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
                x2 = nnz(Plant.Data.Hydro.Timestamp<Date(end) & Plant.Data.Hydro.Timestamp<Plant.Dispatch.Timestamp(1));
                n = x2-x1+1;
                Time(end+1:end+n) = Plant.Data.Hydro.Timestamp(x1:x2); 
                if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                    Flow(end+1:end+n) = Plant.Data.Hydro.OutFlow(x1:x2,I);
                else Flow(end+1:end+n) = Plant.Data.Hydro.SpillFlow(x1:x2,I);
                end
            end
            if any(Date>Plant.Dispatch.Timestamp(1)) 
                x1 = nnz(Plant.Dispatch.Timestamp<Date(1) & Plant.Dispatch.Timestamp>0);
                if x1==0
                    x1 = 1;
                end
                x2 = nnz(Plant.Dispatch.Timestamp<=Date(end) & Plant.Dispatch.Timestamp>0);
                if Plant.Dispatch.Timestamp(x2+1)>0
                    x2 = x2+1;
                end
                Time(end+1:end+x2-x1+1) = Plant.Dispatch.Timestamp(x1:x2);%first time is 1 step before first date, last step is after last date
                Flow(end+1:end+x2-x1+1) = Plant.Dispatch.GeneratorState(x1:x2,nG+nLcum+line);
            end
        end
        riverFlow = interp1(Time,Flow,Date);
    end  
end

