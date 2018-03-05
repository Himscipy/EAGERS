function BuildingWarmUp(Date,days)
global Plant CurrentState
nB = length(Plant.Building);
for i = 1:1:nB
    Build = Plant.Building(i);
    Location = Plant.subNet.Electrical.Location(Build.QPform.Electrical.subnetNode);
    Tzone = 20;
    Twall = 20;
    wuDate = linspace(Date(1),Date(1)+1 - Plant.optimoptions.Resolution/24,24)';
    wuWeather = WeatherForecast(wuDate);
    wuSG = SolarGain(Build,Date,Location,wuWeather);
    wuLoads = BuildingLoads(Build,wuDate,wuSG);
    wuExternalGains = wuSG.Walls + wuSG.Roof;
    for d = 1:1:days
        [~,~,~,Tzone,Twall,~] = BuildingProfile(Build,wuDate,wuLoads.InternalGains,wuExternalGains,wuWeather.Tdb,wuWeather.RH,Tzone(end,1),Twall(end,1));
    end
    if Tzone(end)<10 || Tzone(end)>40
        disp('BuildingWarmUp Error')
    end
    CurrentState.Buildings(1,i) = Tzone(end);
    CurrentState.Buildings(2,i) = Twall(end);
    CurrentState.Buildings(3,i) = Date(1); 
end