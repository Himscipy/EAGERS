function SOC = calculateHydroSOC(Data,ForecastTime)
global Plant CurrentState InOut Last24hour
nS = length(Data(:,1))-1;
nG = length(Plant.Generator);
SOC = zeros(nS+1,nG); %moved from ln 8; was causing to erase SOC
Date = ForecastTime;
dt = (Date(2:end)-Date(1:end-1))*24;

for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        SOC(1,:) = CurrentState.Hydro; 
        DateNow = Plant.Data.HydroHistory.Timestamp(end);
        A = datevec(DateNow);
        hour = ((Date-floor(DateNow))*24);
        n = Plant.Generator(i).QPform.Hydro.subnetNode;
        OutFlow = Data(2:end,nG+Plant.subNet.Hydro.lineNumber(n));
        if strcmp(Plant.optimoptions.forecast,'Perfect')
            %a, what does the historical fit think previous 24 hours of source/sink was
            Forecast = mean(Plant.Data.HistProf.Hydro.SourceSink{i}(A(2),:));
            %b, what was the actual source/sink in that time
            xi = nnz(Plant.Data.HydroHistory.Timestamp<=(DateNow-1));
            Actual = Plant.Data.HydroHistory.SourceSink(xi:end,i); %last 24 hours

            %c, find the scalar ratio of the two
            scale = mean(Actual)/Forecast;

            %d, what does the historical fit think the next ___ hours of source sink will be
            InFlow1(:,i) = interp1(0:24,[Plant.Data.HistProf.Hydro.SourceSink{i}(A(2),end),Plant.Data.HistProf.Hydro.SourceSink{i}(A(2),:)],mod(hour,24));

            %e, scale this expected source sink
            InFlow = InFlow1(:,i)*scale;
%             InFlow(:,1)= InFlow1;
        else  
            InFlow(:,1) = interp1(0:24,[Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),end),Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),:)],mod(hour,24));
        end 
        InOut.SS(:,i) = InFlow;
        K = Plant.subNet.Hydro.UpRiverNodes{n};
        for j = 1:1:length(K)
            T = Plant.subNet.Hydro.lineTime(K(j));   
            OutUpStream = Data(2:end,nG+Plant.subNet.Hydro.lineNumber(K(j)));
            InFlow = InFlow + interp1([Last24hour.Timestamp;Date(2:end)],[Last24hour.Hydro.OutFlow(:,K(j));OutUpStream],[Date-T/24]);
        end 
        
        for t = 1:1:nS
            SOC(t+1,i) = SOC(t,i) + (((sum(InFlow(t))-OutFlow(t)))/12.1).*dt(t);
        end
    end
    InOut.InFlow(:,i) = InFlow(2:end);
    InOut.OutFlow(:,i) = OutFlow;
    InOut.Diff = InOut.InFlow-InOut.OutFlow;

end 

%% Old way to calculate InFlow

%             t = 1;
%             while (Date(t)-DateNow)<T/24
%                 if (Date(t)-DateNow)<T/24
%                     frac = Plant.subNet.Hydro.frac(K(j)); %last time of adding previous outflows, only need a fraction
%                 else frac = 1;
%                 end
%                 InFlow(t) = InFlow(t) + frac*interp1(Plant.Data.HydroHistory.Timestamp,Plant.Data.HydroHistory.OutFlow(:,K(j)),Date(t)-T/24); %outflow from previous dispatches because the time preceeds t = 0
%                 t = t +1;
%             end
%             if frac == 1
%                 InFlow(t) = InFlow(t) + Data(2,nG + Plant.subNet.Hydro.lineNumber(K(j)));
%             else
%             InFlow(t) = InFlow(t) + (1-frac)*Data(2,nG + Plant.subNet.Hydro.lineNumber(K(j)));
%             end
%             tData = 2;
%             while t<nS
%                 t = t+1;
%                 tData = tData+1;
%                 InFlow(t) = InFlow(t) + Data(tData,nG + Plant.subNet.Hydro.lineNumber(K(j)));
%             end
         
%             t0 = hour(1);
%             frac = Plant.subNet.Hydro.frac(K(j)); %last time of adding previous outflows, only need a fraction
%             HistTime = (Plant.Data.HydroHistory.Timestamp(end-ceil(T):end)-Plant.Data.HydroHistory.Timestamp(end))*24; %Last 24 hrs 
%             HistOutFlow = Plant.Data.HydroHistory.OutFlow(:,K(j)); %Last 24 hrs Outflow
%             InFlow = interp1([HistTime,hour(t0+1:end)],[(HistOutFlow(end-ceil(T))*frac),(HistOutFlow(end-floor(T))*(1-frac)),(HistOutFlow(end-floor(T)+1:end)),OutFlow],[HistTime,Date(2:end)]);
        

%% Old way to calc SourceSink

%         x1 = max(1,nnz(Plant.Data.Timestamp<Date(1)));
%             x2 = min(nnz(Plant.Data.Timestamp<Date(end))+1,length(Plant.Data.Timestamp)); 
%             InFlow = interp1(Plant.Data.Timestamp(x1:x2),Plant.Data.Hydro.SourceSink(x1:x2,i),Date);
%        