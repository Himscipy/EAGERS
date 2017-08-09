function Compare2Eplus(building,weather,Date)
% Compare a simulated building to it's energyPlus counterpart
% building is the EAGERS building structure, weather is a sructure of dry
% bulb, wetbuld and relative humididty, Date is the timespan to be
% compared.
%% Load eplus file
global Model_dir
load(fullfile(Model_dir,'Design','Compare2Eplus',strcat(building.Name,'.mat')))           % building

h = gcf;
%% Run EAGERS method
[Equipment,InteriorLighting,ExteriorLighting,Cooling,Heating,Fan_Power,OtherLoads] = BuildingProfile(building,weather,Date);
HVAC_electric = Heating/building.VariableStruct.COP_H + Cooling/building.VariableStruct.COP_C + Fan_Power;
Total = Equipment + InteriorLighting + ExteriorLighting + HVAC_electric + OtherLoads;

D = datevec(Date(1));
DofY = Date - datenum([D(1), 1,1]);
dt = Date(2:end) - Date(1:end-1);

EplusTime = datenum(datevec(mtr.Timestamp));
D = datevec(EplusTime(1));
EplusTime = EplusTime - datenum([D(1), 1,1]);



xi = max(1,nnz((EplusTime)<=DofY(1)));
xf = nnz((EplusTime)<=DofY(end));
Xe = EplusTime(xi:xf);
if xi == 1
    dtE = Xe - ([0;EplusTime(1:xf-1)]);
else dtE = Xe - (EplusTime(xi-1:xf-1));
end

%compare electric (without HVAC)
f1 = figure(1);
set(f1,'name','EquipmentLighting');
h1 = gca;
hold(h1,'off')
E_Plus = mtr.InteriorLightsElectricity(xi:xf) + mtr.InteriorEquipmentElectricity(xi:xf);% + mtr.ExteriorLightsElectricity(xi:xf)
plotStep(Xe,E_Plus)
hold(h1,'on')
plotStep(DofY,Equipment + InteriorLighting,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Equipment(2:end).*dt + InteriorLighting(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in equipment & lighting load is',num2str(PercError)))

%compare cooling electric loads
f2 = figure(2);
set(f2,'name','Cooling');
h2 = gca;
hold(h2,'off')
E_Plus = mtr.CoolingElectricity(xi:xf);
plotStep(Xe,E_Plus)
hold(h2,'on')
plotStep(DofY,Cooling/building.VariableStruct.COP_C,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Cooling(2:end).*dt/building.VariableStruct.COP_C)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in cooling electric load is',num2str(PercError)))

%compare heating electric loads
f3 = figure(3);
set(f3,'name','Heating');
h3 = gca;
hold(h3,'off')
E_Plus = mtr.HeatingElectricity(xi:xf);
plotStep(Xe,E_Plus)
hold(h3,'on')
plotStep(DofY,Heating/building.VariableStruct.COP_H,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Heating(2:end).*dt/building.VariableStruct.COP_H)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in heating electric load is',num2str(PercError)))

%compare net HVAC electric loads
f4 = figure(4);
set(f4,'name','Fans');
h4 = gca;
hold(h4,'off')
E_Plus = mtr.FansElectricity(xi:xf);
plotStep(Xe,E_Plus)
hold(h4,'on')
plotStep(DofY,Fan_Power,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Fan_Power(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in net Fan electric load is',num2str(PercError)))

%compare net HVAC loads
f5 = figure(5);
set(f5,'name','Net HVAC')
h5 = gca;
hold(h5,'off')
E_Plus = mtr.ElectricityHVAC(xi:xf);
plotStep(Xe,E_Plus)
hold(h5,'on')
plotStep(DofY,HVAC_electric,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(HVAC_electric(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in net HVAC electric load is',num2str(PercError)))
%compare exterior lights
E_Plus = mtr.ExteriorLightsElectricity;
PercError = 100*(sum(ExteriorLighting(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in net exterior lighting electric load is',num2str(PercError)))
sum(mtr.ExteriorLightsElectricity(xi:xf))/sum(mtr.ElectricityNetFacility(xi:xf))
%compare net Building electric loads
f6 = figure(6);
set(f6,'name','Net Building')
h6 = gca;
hold(h6,'off')
E_Plus = mtr.ElectricityNetFacility(xi:xf);
plotStep(Xe,E_Plus)
hold(h6,'on')
plotStep(DofY,Total,'r')
xlabel('Day of Year');
ylabel('Electric Load (kW)');
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Total(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
disp(strcat('Total % error in net building electric load is',num2str(PercError)))

figure(h)

end %ends function Compare2ePlus
function plotStep(varargin)
Xe = varargin{1};
Ye = varargin{2};
if length(varargin)>2
    c = varargin{3};
else c = 'b';
end
X = zeros(2*length(Xe),1);
Y = zeros(2*length(Xe),1);
dt = Xe(2:end) - Xe(1:end-1);
dt(end+1) = dt(end);
for t = 1:1:length(Xe)
    X(2*t-1:2*t) = [Xe(t)-dt(t);Xe(t);];
    Y(2*t-1:2*t) = Ye(t);
end
X(end) = Xe(end);
Y(end) = Ye(end);
plot(X,Y,c)
end%Ends function plotEplusStep
