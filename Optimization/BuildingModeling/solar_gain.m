function SG = solar_gain(Build,Date,Location,Weather)
%% Solar gain
nD = length(Date);
or = Build.VariableStruct.Orientation;
sgDiffuse = Build.VariableStruct.WallArea*Weather.irradDiffHorz;%Diffuse irradiance (W)
[SG.Sunrise, SG.Sunset, azimuth, zenith] = solar_calc(Location.Longitude, Location.Latitude, Location.TimeZone, Date);
DN = [cosd(azimuth).*cosd(90 - zenith), sind(azimuth).*cosd(90 - zenith),sind(90 - zenith)];%Direct normal (incoming vector of sunlight)
n = ones(nD,1);
sgDirect = 0.25*Build.VariableStruct.WallArea*Weather.irradDireNorm.*(max(0,dot(n*[cosd(or),sind(or),0],DN,2)) + max(0,dot(n*[-cosd(or),-sind(or),0],DN,2)) + max(0,dot(n*[-sind(or),cosd(or),0],DN,2)) + max(0,dot(n*[sind(or),-cosd(or),0],DN,2)));
SG.Windows = Build.VariableStruct.WindowTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%solar gain (heat) through windows http://www.commercialwindows.org/vt.php
SG.Walls = Build.VariableStruct.WallAbsorption*(1-Build.VariableStruct.WindowWallRatio)*(sgDirect+sgDiffuse)/1000;%solar gain (heat) absorbed by walls
SG.VisibleLight = Build.VariableStruct.LightTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%visible light transmitted through windows
SG.Roof = Build.VariableStruct.RoofArea*Build.VariableStruct.WallAbsorption*Weather.irradDireNorm.*max(0,dot(n*[0,0,1],DN,2))/1000;%solar gain (heat) through roof
end%Ends function solar_gain