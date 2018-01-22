global Plant Last24hour DateSim
Last24hour = [];
%find the current date & load last24hour
if isempty(Plant.Building)
    TimeYesterday = linspace(DateSim-1,DateSim,ceil(24/Plant.optimoptions.Resolution)+1)';
    if Plant.Data.Timestamp(1)<(DateSim-1) && Plant.Data.Timestamp(end)>DateSim
        Last24hour = GetHistoricalData(TimeYesterday);
    else %need to have this in terms of the first timestep
        Last24hour = GetHistoricalData(TimeYesterday+1);
        Last24hour.Timestamp = TimeYesterday;
    end
    Data = GetHistoricalData(DateSim); 
else
    Data = updateForecast(DateSim,[]);
    Date= linspace(DateSim-1,DateSim,24/Plant.optimoptions.Resolution+1)';
    Last24hour = updateForecast(Date,Data);
    %%only doing electric part of building right now. Need to add Cooling
    %%and Heating if Plant has Chillers/Heaters
end