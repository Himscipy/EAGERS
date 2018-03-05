function Power = RenewableOutput(Date,DirNormal)
global Plant
% i is the generator #
%Date is the vector of time points at which you request data
%str chooses between forecast and actual data
nS = length(Date);
nG = length(Plant.Generator);
Power = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Gen = Plant.Generator(i).VariableStruct;
        if isfield(Plant.Generator(i),'QPform')
            Location = Plant.subNet.Electrical.Location(Plant.Generator(i).QPform.Electrical.subnetNode);
        else
            Location = struct('Latitude',40, 'Longitude',-105, 'TimeZone',-7);
        end
        if strcmp(Gen.Type,'Solar') 
            [~,~,Azimuth,Zenith] = SolarCalc(Location.Longitude,Location.Latitude, Location.TimeZone,Date);
            if strcmp(Gen.Tracking,'fixed')
                Power(:,i) = Gen.Sizem2 * (DirNormal/1000).*cosd(Zenith-Gen.Tilt).*(cosd(Azimuth-Gen.Azimuth))*Gen.Eff;
            elseif strcmp(Gen.Tracking,'1axis')
                Power(:,i) = Gen.Sizem2*(DirNormal/1000).*cosd(Zenith-Gen.Tilt)*Gen.Eff;
            else
                Power(:,i) = Gen.Sizem2*(DirNormal/1000)*Gen.Eff; %Dual axis
            end
        elseif strcmp(Gen.Type,'Wind')
            %% need to add wind
        end
    end
end
