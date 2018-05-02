global Model_dir
Model_dir=strrep(which('EAGERS.m'),'EAGERS.m','');
p = path;
dir1 = fullfile(Model_dir,'Projects');
if isempty(strfind(p,dir1))% && isempty(strfind(p,'Projects'))
    SubPaths2Add = {fullfile(Model_dir);
                    fullfile(Model_dir,'Data');
                    fullfile(Model_dir,'Data','Solar');
                    fullfile(Model_dir,'Data','Wind');
                    fullfile(Model_dir,'Data','SavedResults');
                    fullfile(Model_dir,'Data','TestData');
                    fullfile(Model_dir,'Data','TrainedNeuralNetworks');
                    fullfile(Model_dir,'Design');
                    fullfile(Model_dir,'Emissions');
                    fullfile(Model_dir,'GUI');
                    fullfile(Model_dir,'GUI','Design');
                    fullfile(Model_dir,'GUI','Graphics');
                    fullfile(Model_dir,'GUI','Optimization');
                    fullfile(Model_dir,'GUI','Simulation');
                    fullfile(Model_dir,'Model Library');
                    fullfile(Model_dir,'Model Library','Saved Models');
                    fullfile(Model_dir,'Model Library','Saved Linearizations');
                    fullfile(Model_dir,'Optimization');
                    fullfile(Model_dir,'Optimization','AnalysisFunctions');
                    fullfile(Model_dir,'Optimization','BasicFunctions');
                    fullfile(Model_dir,'Optimization','BuildingModeling');
                    fullfile(Model_dir,'Optimization','Communication');
                    fullfile(Model_dir,'Optimization','Forecasting');
                    fullfile(Model_dir,'Optimization','NRELmethod');
                    fullfile(Model_dir,'Optimization','SetupFunctions');
                    fullfile(Model_dir,'Optimization','SolverFunctions');
                    fullfile(Model_dir,'Optimization','RealTimeControl');
                    fullfile(Model_dir,'Optimization','UpdateFunctions');
                    fullfile(Model_dir,'Projects');
                    fullfile(Model_dir,'Simulation');
                    fullfile(Model_dir,'Simulation','BasicFunctions');
                    fullfile(Model_dir,'Simulation','Components');
                    fullfile(Model_dir,'Simulation','CompressorMaps');
                    fullfile(Model_dir,'Simulation','Controls');
                    fullfile(Model_dir,'Simulation','FCMaps');
                    fullfile(Model_dir,'Simulation','LookupFunctions');
                    fullfile(Model_dir,'System Library');
                    fullfile(Model_dir,'System Library','AbChiller');
                    fullfile(Model_dir,'System Library','AC_DC');
                    fullfile(Model_dir,'System Library','AirHeater');
                    fullfile(Model_dir,'System Library','Battery');
                    fullfile(Model_dir,'System Library','Buildings');
                    fullfile(Model_dir,'System Library','Buildings','RealBuildingData');
                    fullfile(Model_dir,'System Library','Chiller');
                    fullfile(Model_dir,'System Library','ColdStor');
                    fullfile(Model_dir,'System Library','Control');
                    fullfile(Model_dir,'System Library','DistrictCool');
                    fullfile(Model_dir,'System Library','DistrictHeat');
                    fullfile(Model_dir,'System Library','Economic');
                    fullfile(Model_dir,'System Library','Electrolyzer');
                    fullfile(Model_dir,'System Library','Grid');
                    fullfile(Model_dir,'System Library','HotStor');
                    fullfile(Model_dir,'System Library','HydrogenStorage');
                    fullfile(Model_dir,'System Library','ICE');
                    fullfile(Model_dir,'System Library','MCFC');
                    fullfile(Model_dir,'System Library','mGT');
                    fullfile(Model_dir,'System Library','NatGas');
                    fullfile(Model_dir,'System Library','PEM');
                    fullfile(Model_dir,'System Library','SOFC');
                    fullfile(Model_dir,'System Library','SolarPV');
                    fullfile(Model_dir,'System Library','SolarStirling');
                    fullfile(Model_dir,'System Library','SolarThermal');
                    fullfile(Model_dir,'System Library','Utility');
                    fullfile(Model_dir,'System Library','WaterHeater');
                    fullfile(Model_dir,'System Library','Weather'); 
                    fullfile(Model_dir,'System Library','Wind');};

    for i=1:length(SubPaths2Add)
        addpath(SubPaths2Add{i});
    end
end
close all
WelcomeScreen