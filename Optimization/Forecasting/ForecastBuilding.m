function [Building,SRtarget] = ForecastBuilding(varargin)
% This function estimates the energy profile of a building 
% Build is a structure of building parameters
% Weather is an hourly weather profile (dry bulb, wet bulb, and relative humidity)
% Date is a vector of points in datenum format at which you are requesting the electric cooling & heating power
global Plant CurrentState
nB = length(Plant.Building);
Weather = varargin{1};
Date = varargin{2};
nS = length(Date);
Building = varargin{3};%forecast with actual loads from TestData
if ~isfield(CurrentState,'Buildings') || abs(round(864000*(CurrentState.Buildings(3,1)+Plant.optimoptions.Resolution/24))/864000 - Date(1))>1e-5
    BuildingWarmUp(Date,6);%%Need warm-up period if not currently running the model
end
Building.E0 = zeros(nS,nB);
Building.C0 = zeros(nS,nB);
Building.H0 = zeros(nS,nB);
Building.Tmin = zeros(nS,nB);
Building.Tmax = zeros(nS,nB);
Building.Tset_H = zeros(nS,nB);
Building.Tset_C = zeros(nS,nB);
Building.Fan_Power = zeros(nS,nB);
Building.AirFlow = zeros(nS,nB);
Building.Tzone = zeros(nS,nB);
Building.Twall = zeros(nS,nB);
Building.Damper = zeros(nS,nB);
for i = 1:1:nB
    Build = Plant.Building(i);
    Tzone = CurrentState.Buildings(1,i);
    Twall = CurrentState.Buildings(2,i);
    [Cooling, Heating, Fan_Power,Tzone,Twall,Damper] = BuildingProfile(Build,Date,Building.InternalGains(:,i),Building.ExternalGains(:,i),Weather.Tdb,Weather.RH,Tzone,Twall);
    Building.Tmin(:,i) = loadSched(Build,Date,'TsetH') - 0.5*Build.VariableStruct.Comfort;
    Building.Tmax(:,i) = loadSched(Build,Date,'TsetC') + 0.5*Build.VariableStruct.Comfort;
    [T_Heat,T_Cool] = EquilibTeperature(Build,Date,Weather,Twall,Building.InternalGains(:,i),Building.ExternalGains(:,i),Building.Tmin(:,i),Building.Tmax(:,i));
    Building.E0(:,i) = Building.NonHVACelectric(:,i) + Fan_Power(:,i);
    Building.C0(:,i) = Cooling(:,i) ;
    Building.H0(:,i) = Heating(:,i);
    if ~Plant.Building(i).QPform.Heating
        %Can upgrade to include non-linear mapping of heating to electricity for built-in HVAC equipment
        Building.E0(:,i) = Building.E0(:,i) + (Building.H0(:,i)- Plant.Building(i).QPform.UA*(Tzone(2:end)-T_Heat))*Build.QPform.H2E;
    end
    if ~Plant.Building(i).QPform.Cooling
        %Can upgrade to include non-linear mapping of cooling to electricity for built-in HVAC equipment
        Building.E0(:,i) = Building.E0(:,i) + (Building.C0(:,i)- Plant.Building(i).QPform.UA*(T_Cool-Tzone(2:end)))*Build.QPform.C2E;
    end
    Building.Tset_H(:,i) = T_Heat;
    Building.Tset_C(:,i) = T_Cool;
    Building.Fan_Power(:,i) = Fan_Power;
    Building.AirFlow(:,i) = Fan_Power/Build.VariableStruct.FanPower;
    Building.Tzone(:,i) = Tzone(2:end);
    Building.Twall(:,i) = Twall(2:end);
    Building.Damper(:,i) = Damper;
end
SRtarget  = 0;
if Plant.optimoptions.SpinReserve
    for i = 1:1:nB
        SRtarget = SRtarget + Building.E0(:,i);
    end
end

end%Ends function ForecastBuilding
