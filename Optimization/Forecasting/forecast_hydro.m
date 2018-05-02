function in_flow = forecast_hydro(prev,date,source_sink,subnet,res)
%forecast the source/sink at each time step: 1 to nS
% does not give a forecast at t=0;
%add the inflow from upstream at t<T
%if there is no connected upstream segment all of the inflow should be in
%the source sink term
n_s = length(date);
in_flow = source_sink;
for n = 1:1:length(subnet.Hydro.nodes) 
    %add the upstream flow
    K = subnet.Hydro.UpRiverNodes{n};
    for j = 1:1:length(K)
        T = subnet.Hydro.lineTime(K(j));
        t = 1;
        while t<=n_s && (date(t)-prev.Timestamp(end))<T/24
            if (date(t)-date(1))< T/24
                frac = subnet.Hydro.frac(K(j)); %last time of adding prevous outflows, only need a fraction
            elseif res>1
                frac = subnet.Hydro.frac(K(j))/(res*24);
            else
                frac = 1;
            end
            in_flow(t,n) = in_flow(t,n) + frac*interp1(prev.Timestamp,prev.Hydro.Outflow(:,K(j)),date(t)-T/24); %outflow from previous dispatches because the time preceeds t = 0
            t = t+1;
        end
    end
end

end %Ends function forecast_hydro
