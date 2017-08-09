function Power = RenewableOutput(i,Date,str)
global Plant DateSim ModelDir
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
    %%Load the Irradiance File, somehow make this a function of location
    if strcmp(str,'Actual') && isfield(Plant.Weather,'Irrad')%% Actual Renewable Generation    
        Irrad = Plant.Weather.Irrad;
    else
        load(fullfile(ModelDir,'Data','DirectNormal.mat'));
        Irrad = DirectNormal(:,1);
    end
    hofY_Irrad = linspace(0,8760,length(Irrad)+1)';
    D = datevec(Date(1));
    h_of_y = mod(24*(Date - datenum([D(1),1,1])),8760); %hour since start of year
    DirNormal = interp1(hofY_Irrad,[Irrad(1);Irrad],h_of_y);
    [~,~,Azimuth,Zenith] = SolarCalc(Location.Longitude,Location.Latitude,Location.TimeZone,Date);
    if strcmp(Gen.Tracking,'fixed')
        Power = Gen.Sizem2*(DirNormal/1000).*cosd(Zenith-Gen.Tilt).*(cosd(Azimuth-Gen.Azimuth))*Gen.Eff;%the data for the renewable power is recorded in 15 min intervals starting at Jan 1st
    elseif strcmp(Gen.Tracking,'1axis')
        Power = Gen.Sizem2*(DirNormal/1000).*cosd(Zenith-Gen.Tilt)*Gen.Eff;
    else Power = Gen.Sizem2*(DirNormal/1000)*Gen.Eff; %Dual axis
    end
elseif strcmp(Gen.Type,'Wind')
    %% need to add wind
end    
