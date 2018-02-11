function ExternalGains = ForecastExternalGains(Forecast)
global Plant
nB = length(Plant.Building);
ExternalGains = zeros(length(Forecast.Timestamp),nB);
%% Solar gain
for i = 1:1:nB
    Build = Plant.Building(i);
    Location = Plant.Building(i).QPform.Location;
    nD = length(Forecast.Timestamp);
    or = Build.VariableStruct.Orientation;
    sgDiffuse = Build.VariableStruct.WallArea*Forecast.Weather.irradDiffHorz;%Diffuse irradiance (W)
    [~, ~, azimuth, zenith] = SolarCalc(Location.Longitude, Location.Latitude, Location.TimeZone,Forecast.Timestamp);
    DN = [cosd(azimuth).*cosd(90 - zenith), sind(azimuth).*cosd(90 - zenith),sind(90 - zenith)];%Direct normal (incoming vector of sunlight)
    n = ones(nD,1);
    sgDirect = 0.25*Build.VariableStruct.WallArea*Forecast.Weather.irradDireNorm.*(max(0,dot(n*[cosd(or),sind(or),0],DN,2)) + max(0,dot(n*[-cosd(or),-sind(or),0],DN,2)) + max(0,dot(n*[-sind(or),cosd(or),0],DN,2)) + max(0,dot(n*[sind(or),-cosd(or),0],DN,2)));
    sgWalls = Build.VariableStruct.WallAbsorption*(1-Build.VariableStruct.WindowWallRatio)*(sgDirect+sgDiffuse)/1000;%solar gain (heat) absorbed by walls
    sgRoof = Build.VariableStruct.RoofArea*Build.VariableStruct.WallAbsorption*Forecast.Weather.irradDireNorm.*max(0,dot(n*[0,0,1],DN,2))/1000;%solar gain (heat) through roof
    ExternalGains(:,i) = sgWalls + sgRoof;
end
end%ends function forecastExternalGains