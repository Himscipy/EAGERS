function [Sunrise,Sunset,Azimuth,Zenith] = SolarCalc(Long,Lat,TimeZone,Date)
%[Sunrise,Sunset] = SunriseSunset(Long, Lat,TimeZone,Date)
% calculated using NOAA solar calculations available at:
% https://www.esrl.noaa.gov/gmd/grad/solcalc/NOAA_Solar_Calculations_day.xls
% Long is longitude where + is to the east
% Lat is latitude with + to the north pole
% TimeZone is the time zone with + to the east
% Date is the date in matlab datenum format, ie Jan 1 2017 = 736696;
% Sunrise and sunset are given in fraction of the day, ie 6am = 6/24
% Azimuth and Zenith are given in degrees
TPM = Date - floor(Date);
JD = Date + 1721058.5 + TPM - TimeZone/24;
JC = (JD-2451545)/36525;
GeomMeanLongSun = mod(280.46646+JC.*(36000.76983+JC*.0003032),360);%degrees
GeomMeanAnom = 357.52911+JC.*(35999.05029 - 0.0001537*JC);%degrees
EccentEarthOrbit = 0.016708634-JC.*(0.000042037+0.0000001267*JC);
SunEqOfCentr = sin(GeomMeanAnom/57.2958).*(1.914602-JC.*(0.004817+0.000014*JC))+sin(2*GeomMeanAnom/57.2958).*(0.019993-0.000101*JC)+sin(3*GeomMeanAnom/57.2958)*0.000289;
SunTrueLong = GeomMeanLongSun + SunEqOfCentr;%degrees
MeanObliqueCorr = 23+(26+((21.448-JC.*(46.815+JC.*(.00059-JC*0.001813))))/60)/60; %degrees
SunAppLong = SunTrueLong - 0.00569-0.00478*sin((125.04-1934.136*JC)/57.2958); %degrees
ObliqueCorr = MeanObliqueCorr + 0.00256*cos((125.04-1934.136*JC)/57.2958); %degrees
declination = asin(sin(ObliqueCorr/57.2958).*sin(SunAppLong/57.2958));%Solar declination in radians
var_y = tan(ObliqueCorr/2/57.2958).*tan(ObliqueCorr/2/57.2958);
EqOfTime = 4*360/(2*pi)*(var_y.*sin(2*GeomMeanLongSun/57.2958)-2*EccentEarthOrbit.*sin(GeomMeanAnom/57.2958)+4*EccentEarthOrbit.*var_y.*sin(GeomMeanAnom/57.2958).*cos(2*GeomMeanLongSun/57.2958)-.5*var_y.*var_y.*sin(4*GeomMeanLongSun/57.2958)-1.25*EccentEarthOrbit.*EccentEarthOrbit.*sin(2*GeomMeanAnom/57.2958));
HAsunrise = 360/(2*pi)*acos(cos(90.833/57.2958)./(cos(Lat/57.2958).*cos(declination))-tan(Lat/57.2958).*tan(declination));%Degrees (sunlight hours)
SolarNoon = (720-4*Long - EqOfTime + TimeZone*60)/1440;
Sunrise = SolarNoon-HAsunrise*4/1440;
Sunset = SolarNoon+HAsunrise*4/1440;

TST = mod(TPM*1440+EqOfTime+4*Long-60*TimeZone,1440); %TrueSolar Time in min
if (TST/4)<0
    HourAngle = TST/4+180; %degrees
else HourAngle = TST/4-180; %degrees
end
Zenith = 360/(2*pi)*acos(sin(Lat/57.2958).*sin(declination)+cos(Lat/57.2958).*cos(declination).*cos(HourAngle/57.2958));%Degrees
if HourAngle >0
    Azimuth = mod(360/(2*pi)*acos(((sin(Lat/57.2958).*cos(Zenith/57.2958))-sin(declination))./(cos(Lat/57.2958).*sin(Zenith/57.2958)))+180,360);%Degrees clockwise from N
else
    Azimuth = mod(540-360/(2*pi)*real(acos(((sin(Lat/57.2958).*cos(Zenith/57.2958))-sin(declination))./(cos(Lat/57.2958).*sin(Zenith/57.2958)))),360);%Degrees clockwise from N
end
end %Ends function SunriseSunset

