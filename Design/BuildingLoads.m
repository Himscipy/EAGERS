function [InternalGains,ExternalGains,Equipment, InteriorLighting, ExteriorLighting, OtherLoads] = BuildingLoads(Build,irradDireNorm,irradDiffHorz,Location,Date)
ExtLights_B = .04;%.125; % fraction of day exterior lights remain on after sunrise
% Time
[Y,Month,D,H,M,S] = datevec(Date);
h1 = H(1) + M(1) + S(1);
D(H==0&M==0&S==0) = D(H==0&M==0&S==0) - 1;
if h1 == 0 % put start date back.
    D(1) = D(1) + 1;
end
days = ceil(Date(end) - datenum([Y(1),Month(1),D(1)]));
sd = floor(Date(1)); % start date
nS = length(Date); % number timesteps

%% Fill out schedule values for each timestep in input Date
% If next point in schedule is shorter than ramp, then ramp shortens to half the gap
p = 1;
DateIndex = cell(days);
for d = 1:1:days
    p2 = p;
    while p2<nS && Date(p2+1)<=sd+d
        p2 = p2+1;
    end
    DateIndex(d) = {p:p2};
    p = p2+1;
end
% Schedules
sched = fieldnames(Build.Schedule);
for i = 1:1:length(sched)
    Profile.(sched{i}) = loadSched(Build,Date,sched{i}); %interp1(Out.(sched{i})(:,1),Out.(sched{i})(:,2), H(p:p2)+M(p:p2)/60+S(p:p2)/3600);
end
%% Equipment load values
Equipment = Build.Area*Build.VariableStruct.equipment/1000*Profile.equipment; % kW of equipment load

%% Other load values
OtherLoads = zeros(nS,1);
if isfield(Build.VariableStruct,'WaterSystem')
    OtherLoads = Build.Area*Build.VariableStruct.WaterSystem*Profile.watersystem;
end

%% Exterior lighting load values
if isfield(Profile,'exteriorlights_solarcontroled')
    for day = 1:days
        currDate = floor(Date(DateIndex{day}(1)));
        HoD = 24*(Date(DateIndex{day}) - currDate); % hour of day
        frac_dark = ones(length(HoD),1);
        [Sunrise, Sunset, ~, ~] = SolarCalc(Location.Longitude,Location.Latitude, Location.TimeZone, currDate);
        
        % Find hour at which sunrise occurs
        h_sr = 1;
        while HoD(h_sr)<(Sunrise*24+ExtLights_B) && h_sr<length(HoD)
            h_sr = h_sr+1;
        end
        if h_sr == 1
            frac_dark(1) = (Sunrise*24+ExtLights_B)/HoD(1);
        else
            frac_dark(h_sr) = ((Sunrise*24+ExtLights_B) - HoD(h_sr-1))/(HoD(h_sr)-HoD(h_sr-1));
        end
        
        % Find hour at which sunset occurs
        h_ss = h_sr+1;
        if h_ss < length(HoD)
            while HoD(h_ss)<(Sunset*24-ExtLights_B)  && h_ss<length(HoD)
                h_ss = h_ss+1;
            end
            frac_dark(h_ss) = 1-((Sunset*24-ExtLights_B) - HoD(h_ss-1))/(HoD(h_ss)-HoD(h_ss-1));
        end
        frac_dark(h_sr+1:h_ss-1) = 0;
        Profile.exteriorlights_solarcontroled(DateIndex{day}) = min(Profile.exteriorlights_solarcontroled(DateIndex{day}'),frac_dark);
    end
    ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights_fixed + Build.VariableStruct.ExteriorLights*Profile.exteriorlights_solarcontroled;
else
    ExteriorLighting = Build.VariableStruct.ExteriorLights*Profile.exteriorlights;
end


%% Solar gain
nD = length(Date);
or = Build.VariableStruct.Orientation;
sgDiffuse = Build.VariableStruct.WallArea*irradDiffHorz;%Diffuse irradiance (W)
[~, ~, azimuth, zenith] = SolarCalc(Location.Longitude, Location.Latitude, Location.TimeZone, Date);
DN = [cosd(azimuth).*cosd(90 - zenith), sind(azimuth).*cosd(90 - zenith),sind(90 - zenith)];%Direct normal (incoming vector of sunlight)
n = ones(nD,1);
sgDirect = 0.25*Build.VariableStruct.WallArea*irradDireNorm.*(max(0,dot(n*[cosd(or),sind(or),0],DN,2)) + max(0,dot(n*[-cosd(or),-sind(or),0],DN,2)) + max(0,dot(n*[-sind(or),cosd(or),0],DN,2)) + max(0,dot(n*[sind(or),-cosd(or),0],DN,2)));
sgWindows = Build.VariableStruct.WindowTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%solar gain (heat) through windows http://www.commercialwindows.org/vt.php
sgWalls = Build.VariableStruct.WallAbsorption*(1-Build.VariableStruct.WindowWallRatio)*(sgDirect+sgDiffuse)/1000;%solar gain (heat) absorbed by walls
vlWindows = Build.VariableStruct.LightTransmittance*Build.VariableStruct.WindowWallRatio*(sgDirect+sgDiffuse)/1000;%visible light transmitted through windows
sgRoof = Build.VariableStruct.RoofArea*Build.VariableStruct.WallAbsorption*irradDireNorm.*max(0,dot(n*[0,0,1],DN,2))/1000;%solar gain (heat) through roof

%% Interior lighting load values
InteriorLighting = Build.Area*Build.VariableStruct.InteriorLights/1000*Profile.interiorlights.*(1 - Build.VariableStruct.DaylightPercent*vlWindows/max(vlWindows)); % kW of lighting load

%% Internal gains
internalGainOccupants = Build.Area * Build.VariableStruct.occupancy * Profile.occupancy * 0.120; % heat from occupants (120 W)
InternalGains = internalGainOccupants + Equipment + InteriorLighting + sgWindows;
ExternalGains = sgWalls + sgRoof;

if isfield(Build.VariableStruct,'DataCenter')%data center is an equipment load but not an internal gain into the zones
    Equipment = Equipment + Build.VariableStruct.DataCenter*Profile.datacenter;
end
end%Ends function Building Loads
