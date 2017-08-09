function InFlow = ForecastHydro(Date,t0)
global Plant
%forecast the source/sink at each time step: 1 to nS
% does not give a forecast at t=0;
%add the inflow from upstream at t<T
%if there is no connected upstream segment all of the inflow should be in
%the source sink term
A = datevec(t0);
hour = (Date-floor(t0))*24;
nS = length(Date);
InFlow = zeros(nS,length(Plant.subNet.Hydro.nodes));
for n = 1:1:length(Plant.subNet.Hydro.nodes) 
    %a, what does the historical fit think previous 24 hours of source/sink was
    Forecast = mean(Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),:));
    %b, what was the actual source/sink in that time
    xi = nnz(Plant.Data.HydroHistory.Timestamp<=(t0-1));
    Actual = Plant.Data.HydroHistory.SourceSink(xi:end,n); %last 24 hours

    %c, find the scalar ratio of the two
    scale = mean(Actual)/Forecast;

    %d, what does the historical fit think the next ___ hours of source sink will be
    InFlow(:,n) = interp1(0:24,[Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),end),Plant.Data.HistProf.Hydro.SourceSink{n}(A(2),:)],mod(hour,24));

    %e, scale this expected source sink
    InFlow(:,n) = InFlow(:,n)*scale;

    %f, add the upstream flow
    K = Plant.subNet.Hydro.UpRiverNodes{n};
    for j = 1:1:length(K)
        T = Plant.subNet.Hydro.lineTime(K(j));
        t = 1;
        while t<=nS &&(Date(t)-t0)<T/24
            if (Date(t)-t0)<T/24
                frac = Plant.subNet.Hydro.frac(K(j)); %last time of adding prevous outflows, only need a fraction
            else frac = 1;
            end
            InFlow(t,n) = InFlow(t,n) + frac*interp1(Plant.Data.HydroHistory.Timestamp,Plant.Data.HydroHistory.OutFlow(:,K(j)),Date(t)-T/24); %outflow from previous dispatches because the time preceeds t = 0
            t = t+1;
        end
    end
end
end%Ends function ForecastHydro
