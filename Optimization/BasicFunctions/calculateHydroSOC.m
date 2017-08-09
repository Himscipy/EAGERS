function SOC = calculateHydroSOC(Data)
global Plant CurrentState
nS = length(Data(:,1))-1;
nG = length(Plant.Generator);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        SOC = zeros(nS+1,nG);
        SOC(1,:) = CurrentState.Hydro; 
        DateNow = Plant.Data.HydroHistory.Timestamp(end);
        A = datevec(DateNow);
        hour = (Date-floor(DateNow))*24;
        n = Plant.Generator(i).QPform.subnetNode;
        OutFlow = Data(2:end,nG+Plant.subNet.Hydro.LineNumber(n));
        InFlow = interp1(0:24,[Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),end),Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),:)],mod(hour,24));
        K = Plant.subNet.Hydro.UpRiverNodes{n};
        for j = 1:1:length(K)
            T = Plant.subNet.Hydro.lineTime(K(j));
            t = 1;
            while (Date(t)-DateNow)<T/24
                if (Date(t)-DateNow)<T/24
                    frac = Plant.subNet.Hydro.frac(K(j)); %last time of adding prevous outflows, only need a fraction
                else frac = 1;
                end
                InFlow(t,n) = InFlow(t,n) + frac*interp1(Plant.Data.HydroHistory.Timestamp,Plant.Data.HydroHistory.OutFlow(:,K(j)),Date(t)-T/24); %outflow from previous dispatches because the time preceeds t = 0
            end
            InFlow(t,n) = InFlow(t,n) + (1-frac)*Data(2,nG + K(j));
            tData = 2;
            while t<nS
                t = t+1;
                tData = tData+1;
                InFlow(t,n) = InFlow(t,n) + Data(tData,nG + K(j));
            end
        end
        for t = 1:1:nS
            SOC(t+1,i) = SOC(t,i) + (sum(InFlow(t))-OutFlow(t))/12.1;
        end
    end
end