function InFlow = ForecastHydro(Date,SourceSink)
global Plant
%forecast the source/sink at each time step: 1 to nS
% does not give a forecast at t=0;
%add the inflow from upstream at t<T
%if there is no connected upstream segment all of the inflow should be in
%the source sink term
global Last24hour
nS = length(Date);
InFlow = SourceSink;
for n = 1:1:length(Plant.subNet.Hydro.nodes) 
    %add the upstream flow
    K = Plant.subNet.Hydro.UpRiverNodes{n};
    for j = 1:1:length(K)
        T = Plant.subNet.Hydro.lineTime(K(j));
        t = 1;
        while t<=nS && (Date(t)-Last24hour.Timestamp(end))<T/24
            if (Date(t)-Date(1))< T/24
                frac = Plant.subNet.Hydro.frac(K(j)); %last time of adding prevous outflows, only need a fraction
            elseif Plant.optimoptions.Resolution/24>1
                frac = Plant.subNet.Hydro.frac(K(j))/Plant.optimoptions.Resolution;
            else
                frac = 1;
            end
            InFlow(t,n) = InFlow(t,n) + frac*interp1(Last24hour.Timestamp,Last24hour.Hydro.OutFlow(:,K(j)),Date(t)-T/24); %outflow from previous dispatches because the time preceeds t = 0
            t = t+1;
        end
    end
end

end %Ends function ForecastHydro
