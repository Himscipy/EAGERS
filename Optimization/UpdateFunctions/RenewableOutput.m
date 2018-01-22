function Power = RenewableOutput(i,Date,str,DirNormal)
global Plant DateSim Model_dir
% i is the generator #
%Date is the vector of time points at which you request data
%str chooses between forecast and actual data
Gen = Plant.Generator(i).VariableStruct;
Enode = Plant.Generator(i).QPform.Electrical.subnetNode;
Location = Plant.subNet.Electrical.Location(Enode);
if isempty(Date)
    Date = DateSim;
end
if strcmp(Gen.Type,'Solar') 
    if isempty(DirNormal)
        %% Load the Irradiance File, somehow make this a function of location
        if isfield(Plant.Weather,'irradDireNorm') 
            Irrad = Plant.Weather.irradDireNorm;
        else
            load(fullfile(Model_dir,'Data','Solar','DirectNormal.mat'));
            Irrad = DirectNormal(:,1);
        end
        hofY_Irrad = linspace(0,8760,length(Irrad)+1)';
        D = datevec(Date(1));
        if strcmp(str,'Actual') || strcmp(Plant.optimoptions.method,'Planning') %make entire forecast perfect if planning
            h_of_y = mod(24*(Date - datenum([D(1),1,1])),8760); % Hour since start of year
        else %use yesterday's actual
            h_of_y = mod(24*(Date - 1 - datenum([D(1),1,1])),8760); % Hour since start of year (yesterday)
        end
        DirNormal = interp1(hofY_Irrad,[Irrad(1);Irrad],h_of_y);
    end
    [~,~,Azimuth,Zenith] = SolarCalc(Location.Longitude,Location.Latitude, Location.TimeZone,Date);
    if strcmp(Gen.Tracking,'fixed')
        Power = Gen.Sizem2 * (DirNormal/1000).*cosd(Zenith-Gen.Tilt).*(cosd(Azimuth-Gen.Azimuth))*Gen.Eff;
    elseif strcmp(Gen.Tracking,'1axis')
        Power = Gen.Sizem2*(DirNormal/1000).*cosd(Zenith-Gen.Tilt)*Gen.Eff;
    else
        Power = Gen.Sizem2*(DirNormal/1000)*Gen.Eff; %Dual axis
    end
elseif strcmp(Gen.Type,'Wind')
    %% need to add wind
end
