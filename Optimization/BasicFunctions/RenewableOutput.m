function Ren= RenewableOutput(i,Date,str)
global Plant DateSim
% i is the generator #
%Date is the vector of time points at which you request data
%str chooses between forecast and actual data
Gen = Plant.Generator(i).VariableStruct;
if isempty(Date)
    Date = DateSim;
end
A = datevec(Date(1));
B = datenum(A(1),1,1);
day1 = floor(Date(1)-B);
day2 = ceil(Date(end)-B);
dataDate = (day1:1/96:day2)+B;
if strcmp(str,'Predict')
    day1 = day1-4;
    day2 = day2+3;    
end
if day1 ==day2 %only happens for a single point, right at midnight
    dataDate = [Date Date+1/96];
    index = [max(1,day1*96),day1*96+1];
elseif day1<=0
    index = [(day1+365)*96:365*96, 1:day2*96];
elseif day2<=365
    index = [day1*96:day2*96];
else
    index = [day1*96:365*96, 1:(day2-365)*96];
end
if strcmp(Gen.Type,'Solar') 
    if strcmp(Gen.Tracking,'fixed')
        Power = Gen.Sizem2*(Gen.Irrad(index)/1000).*cosd(Gen.SunZen(index)-Gen.Tilt).*(cosd(Gen.SunAz(index)-Gen.Azimuth))*Gen.Eff;%the data for the renewable power is recorded in 15 min intervals starting at Jan 1st
    elseif strcmp(Gen.Tracking,'1axis')
        Power = Gen.Sizem2*(Gen.Irrad(index)/1000).*cosd(Gen.SunZen(index)-Gen.Tilt)*Gen.Eff;
    else Power = Gen.Sizem2*(Gen.Irrad(index)/1000)*Gen.Eff; %Dual axis
    end
elseif strcmp(Gen.Type,'Wind')
    %% need to add wind
end    
if strcmp(str,'Actual')%% Actual Renewable Generation     
    Ren = interp1(dataDate,Power,Date);
elseif strcmp(str,'Predict')
    AveragedDays = (Power(1:end-7*96)+Power(97:end-6*96) + Power(2*96+1:end-5*96) + Power(3*96+1:end-4*96) + Power(4*96+1:end-3*96) + Power(5*96+1:end-2*96) + Power(6*96+1:end-1*96) + Power(7*96+1:end))/8; %average 4 prior days and next 3 days
    Ren = interp1(dataDate, AveragedDays, Date);
end
if any(isnan(Ren))
    disp('WTF')
end