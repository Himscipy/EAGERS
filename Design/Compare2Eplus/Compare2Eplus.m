function errors = Compare2Eplus(building,weather,Date)
% Compare a simulated building to it's energyPlus counterpart
% building is the EAGERS building structure, weather is a sructure of dry
% bulb, wetbuld and relative humididty, Date is the timespan to be
% compared.

%% Constants
N_FIGS = 8;

%% Handle input DATE
[s1, ~] = size(Date);
if s1 == 1
    Date = Date';
end

%% Load EPlus meter file
global Model_dir
load(fullfile(Model_dir,'Design','Compare2Eplus',strcat(building.Name,'.mat')))           % building

%% Run EAGERS method
[Equipment, InteriorLighting, ExteriorLighting, Cooling, Heating, Fan_Power, ...
    OtherLoads] = BuildingProfile(building, weather, Date);
HVAC_electric = Heating/building.VariableStruct.COP_H + ...
    Cooling/building.VariableStruct.COP_C + Fan_Power;
Total = Equipment + InteriorLighting + ExteriorLighting + HVAC_electric + ...
    OtherLoads;

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
else
    dtE = Xe - (EplusTime(xi-1:xf-1));
end

%% Plotting
errors = zeros(N_FIGS, 1);
fprintf('Total %% errors:\n')

%compare equipment
figTitle = 'Equipment';
figNum = 1;
ePlus = mtr.InteriorEquipmentElectricity(xi:xf);
eagers = Equipment;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare interior lighting
figTitle = 'Interior Lighting';
figNum = 2;
ePlus = mtr.InteriorLightsElectricity(xi:xf);
eagers = InteriorLighting;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare exterior lighting
figTitle = 'Exterior Lighting';
figNum = 3;
ePlus = mtr.ExteriorLightsElectricity(xi:xf);
eagers = ExteriorLighting;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare cooling electric loads
figTitle = 'Cooling';
figNum = 4;
ePlus = mtr.CoolingElectricity(xi:xf);
eagers = Cooling / building.VariableStruct.COP_C;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare heating electric loads
figTitle = 'Heating';
figNum = 5;
ePlus = mtr.HeatingElectricity(xi:xf);
eagers = Heating / building.VariableStruct.COP_H;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare fan loads
figTitle = 'Fans';
figNum = 6;
ePlus = mtr.FansElectricity(xi:xf);
eagers = Fan_Power;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare net HVAC loads
figTitle = 'Net HVAC';
figNum = 7;
ePlus = mtr.ElectricityHVAC(xi:xf);
eagers = HVAC_electric;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

%compare net Building electric loads
figTitle = 'Net Building';
figNum = 8;
ePlus = mtr.ElectricityNetFacility(xi:xf);
eagers = Total;
errors(figNum) = ...
    elec_load_plot_setup(figTitle, figNum, Xe, ePlus, DofY, eagers, dt, dtE);

for i = N_FIGS:-1:1
    figure(i)
end

end % Ends function Compare2Eplus


function PercError = elec_load_plot_setup(figTitle, figNum, Xe, E_Plus, ...
    DofY, Eagers, dt, dtE)

fig = figure(figNum);
set(fig,'name',figTitle)
hold off
plotStep(Xe,E_Plus,'b')
hold on
plotStep(DofY,Eagers,'r')
xlabel('Day of Year')
ylabel('Electric Load (kW)')
legend({'EnergyPlus','EAGERS'})
PercError = 100*(sum(Eagers(2:end).*dt)-sum(E_Plus.*dtE))/sum(E_Plus.*dtE);
fprintf('%s electric load:\t%f\n', figTitle, PercError)

end % Ends function elec_load_plot_setup


function plotStep(varargin)
Xe = varargin{1};
Ye = varargin{2};
if length(varargin) > 2
    lineSpec = {varargin{3:end}};
else
    lineSpec = {};
end

X = [0; Xe];
Y = [Ye; Ye(end)];
stairs(X, Y, lineSpec{:})
end % Ends function plotStep
