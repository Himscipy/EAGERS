function buildings = building_warmup(buildings,wu_weather,res,date,days)
n_b = length(buildings);
for i = 1:1:n_b
    building_i = buildings(i);
    zone = 20;
    wall = 20;
    if length(date) == 1
        wuDate = linspace(date(1) + res,date(1)+1,24)';
    else
        wuDate = linspace(date(1),date(1)+1 - res,24)';
    end
    wu_sg = solar_gain(building_i,wuDate,building_i.QPform.Location,wu_weather);
    wu_loads = building_loads(building_i,wuDate,wu_sg);
    wu_external_gains = wu_sg.Walls + wu_sg.Roof;
    for d = 1:1:days
        [~,~,~,zone,wall,~] = building_profile(building_i,wuDate,wu_loads.InternalGains,wu_external_gains,wu_weather.Tdb,wu_weather.RH,zone(end,1),wall(end,1));
    end
    if zone(end)<10 || zone(end)>40
        disp('BuildingWarmUp Error')
    end
    buildings(i).Tzone = zone(end);
    buildings(i).Twall = wall(end);
    buildings(i).Timestamp = date(1);
end
end%Ends function building_warmup