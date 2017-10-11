function solarGain = SolarGainThroughWindows(building, location, weather, date)
%SOLARGAINTHROUGHWINDOWS Calculates solar gain through a building's
%windows.
%   solarGain = SolarGainThroughWindows(building,weather,date)
%   Returns vector of solar gain values corresponding to each time in
%   datetime array DATE for weather conditions specified by struct WEATHER.
%   BUILDING is a struct containing information about the building for
%   which the solar gain should be calculated. The solar gain is calculated
%   using the dot product of window area and solar irradiance vectors. The
%   return value units will correspond to the units of IRRADIANCE
%   multiplied by m^2. (E.g. IRRADIANCE in Wh/m^2 will result in return
%   value in Wh.)
%
%   ASSUMPTIONS:
%   --  Building footprint is square, with four exterior walls of equal
%       dimensions, and window area is evenly distributed among those four
%       walls.
%   --  Coordinate system is Cartesian, with origin at southwestern-most
%       corner of the building. x corresponds with north, y corresponds
%       with east, and the z-axis is the zenith.
%   --  All angles in degrees.
%   --  The input DATE contains irradiance values of zero at times when the
%       sun has not yet risen.
%   --  The input DATE does not include any values for Feb 29.
%   --  Position of sun is calculated by function SolarCalc.
%
%   See also SOLARCALC

wallAreaGross   = building.VariableStruct.ExteriorWallAreaGross;
winWallRatio    = building.VariableStruct.WindowWallRatio;
totalWinArea    = wallAreaGross * winWallRatio;
orientation     = building.VariableStruct.Orientation;

areaVectors     = window_area_vectors(totalWinArea, orientation);
sgDirect        = direct_solar_gain(areaVectors, weather.irradDireNorm, ...
                    location, date);
sgDiffuse       = diffuse_solar_gain(totalWinArea, weather.irradDiffHorz);

solarGain = sgDirect + sgDiffuse;

end


%% Internal functions

% Area vectors for all four cardinal directions.
function area = window_area_vectors(totalWinArea,orientation)

orCos = cosd(orientation);
orSin = sind(orientation);

aNorth  = 1/4 * totalWinArea * [orCos   orSin   0];
aSouth  = -aNorth;
aEast   = 1/4 * totalWinArea * [-orSin  orCos   0];
aWest   = -aEast;

area = [aNorth; aSouth; aEast; aWest];

end


% Solar gains for all four cardinal directions.
function sgDirect = direct_solar_gain(areaVecs, vecDireNorm, loc, vecDate)

sizeDate = size(vecDate);
sgDirect = zeros(sizeDate);
for i = 1:1:length(vecDate)
    sgDirect(i)  = solar_gain(areaVecs, vecDireNorm(i), loc, vecDate(i));
end

end


% Solar gains for all four cardinal directions.
function sgDiffuse = diffuse_solar_gain(totWinArea, vecDiffHorz)

sizeDiff = size(vecDiffHorz);
sgDiffuse = zeros(sizeDiff);
for i = 1:1:length(vecDiffHorz)
    sgDiffuse(i)  = vecDiffHorz(i) * totWinArea;
end

end


% Solar gain for a single area vector, irradiance value, location and date.
function sGain = solar_gain(areaVecs, irrad, loc, date)

condition = size(date)==[1,1];
assert(sum(condition)==2, 'Input DATE must be a single value.')

lat         = loc.Latitude;
long        = loc.Longitude;
timeZone    = loc.TimeZone;
[~, ~, azimuth, zenith] = SolarCalc(long, lat, timeZone, date);

% Irradiance vector
altitude = 90 - zenith;
xSun = irrad * cosd(azimuth) * cosd(altitude);
ySun = irrad * sind(azimuth) * cosd(altitude);
zSun = irrad * sind(altitude);
directNormalIrradiance = [xSun, ySun, zSun];

% Direct irradiance only happens between sunrise and sunset
sGain = 0;
if sum(abs(directNormalIrradiance)) > 0
    nDirs = 4; % N, S, E, W
    for i = 1:1:nDirs
        sgDirection = dot(areaVecs(i,:), directNormalIrradiance);
        % Solar gain is never negative
        sGain = sGain + max(sgDirection, 0);
    end
end

end
