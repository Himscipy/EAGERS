% This function loads the test data that the projects will be run against
%Option a) use simulated building if anything in Plant.Building (need to double check for other needed data streams like solar/wind/hydro)
%Option b) use historical record already in Plant.Data
%Option c) load test data from a file
global Plant TestData Model_dir
if isempty(TestData)
    %check for Plant.Building
    if isfield(Plant,'Data') && ~isempty(Plant.Data)
        K = center_menu('Select Test Data', 'Simulated Building', 'Historical Data in project', 'Load from File');
    else
        K = center_menu('Select Test Data', 'Simulated Building', 'Load from File');
        if K ==2
            K = 3;
        end
    end
    if K == 1
        if isfield(Plant,'Building') && ~isempty(Plant.Building)
            K2 = center_menu('Select Building', 'Building(s) in Project', 'Load from File');
            if K2 == 2
                files = dir(fullfile(Model_dir, 'System Library','Buildings','*.mat'));
                listB=strrep({files.name},'.mat','');
                [s,v] = listdlg('ListString',listB,'SelectionMode','multiple','Name','Saved buildings','PromptString','Select one or more buildings');
                Plant.Building = [];
                for i = 1:1:length(s)
                    load(fullfile(Model_dir, 'System Library','Buildings',strcat(listB{s(i)},'*.mat')));
                    Plant.Building(i) = building;
                end
            end
        end
        %% select simulation dates
        nS = 365*24/Plant.optimoptions.Resolution+1;
        TestData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),nS)';
        %load weather file if not already loaded
        cd(fullfile(Model_dir,'System Library','Weather'))
        [fn,pn,~] = uigetfile('*.mat','Load Weather File');
        load(fullfile(pn,fn));
        cd(Model_dir)
        TestData.Weather = interpolateWeather(weather,TestData.Timestamp);
        reLoadBuilding(Plant.Building,Plant.Network)
    elseif K == 2
        TestData = Plant.Data;%% Revise this so you can pull from more than what is loaded in Plant
        if isfield(TestData,'HistProf')
            TestData = rmfield(TestData,'HistProf');
        end
        if isfield(Plant.Data,'Hydro')
            TestData = rmfield(TestData,'Hydro');
            TestData.Hydro.SourceSink = Plant.Data.Hydro.SourceSink;
            TestData.Hydro.OutFlow = Plant.Data.Hydro.OutFlow;
        end  
    elseif K == 3
        cd(fullfile(Model_dir,'Data','TestData'))
        [fn,pn,~] = uigetfile('*.mat','Load Test Data File');
        load(fullfile(pn,fn));
        cd(Model_dir)
    end
end