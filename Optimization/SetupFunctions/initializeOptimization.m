global Plant Model_dir TestData
if ~isfield(TestData,'Weather') || ~isfield(TestData.Weather,'irradDireNorm') %no solar data
    nG = length(Plant.Generator);
    solar = false;
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'Solar')
            solar = true;
        end
    end
    if solar %need to load solar profile
        cd(fullfile(Model_dir,'Data','Solar'))
        [fn,pn,~] = uigetfile('*.mat','Load Solar Irradiation Data File');
        load(fullfile(pn,fn));
        TestData.Weather.irradDireNorm = zeros(length(TestData.Timestamp),1);
        if TestData.Timestamp(1)>=weather.Timestamp(1) && TestData.Timestamp(end)<=weather.Timestamp(end)
            TestData.Weather.irradDireNorm = interp1(weather.Timestamp,weather.irradDireNorm,TestData.Timestamp); 
        else %need to move solar data timestamp to line up with TestData timestamp
            D1 = datevec(TestData.Timestamp(1));
            D2 = datevec(weather.Timestamp(1));
            if datenum([D1(1),D2(2),D2(3),D2(4),D2(5),D2(6)])>TestData.Timestamp(1)
                D2(1) = D2(1)-1;
            end
            weather.Timestamp = weather.Timestamp + datenum([D1(1),1,1]) - datenum([D2(1),1,1]);%days between start of years for testdata and solar data respectively
            y = 1;
            Xi = 1;
            Xf = nnz(TestData.Timestamp<=weather.Timestamp(end));
            while Xi<=length(TestData.Timestamp)
                TestData.Weather.irradDireNorm(Xi:Xf,1) = interp1(weather.Timestamp,weather.irradDireNorm,TestData.Timestamp(Xi:Xf)); 
                weather.Timestamp = weather.Timestamp + datenum([D1(1)+y,1,1]) - datenum([D1(1)+y-1,1,1]);%shift weather data 1 year
                Xi = Xf+1;
                Xf = nnz(TestData.Timestamp<=weather.Timestamp(end));
                y = y+1;
            end
        end
        cd(Model_dir)
    end
end

if ~isfield(Plant,'Data')
    if isfield(TestData,'Weather')
        S = fieldnames(TestData.Weather);
        for j = 1:1:length(S)
            if isnumeric(TestData.Weather.(S{j}))
                Plant.Data.HistProf.(S{j}) = TypicalDay([],TestData.Timestamp,TestData.Weather.(S{j}));
            end
        end
    end
    if isfield(TestData,'Hydro')
        nodes = length(TestData.Hydro.SourceSink(1,:));
        for n = 1:1:nodes
            Plant.Data.HistProf.Hydro.SourceSink(n) = {TypicalDay([],TestData.Hydro.Timestamp,TestData.Hydro.SourceSink(:,n))};
        end
    end
elseif ~isfield(Plant.Data,'HistProf') || isempty(Plant.Data.HistProf)
    S = fieldnames(TestData.Weather);
    for j = 1:1:length(S)
        if isnumeric(TestData.Weather.(S{j}))
            if isfield(Plant.Data,'Weather') && isfield(Plant.Data.Weather,S{j})
                Plant.Data.HistProf.(S{j}) = TypicalDay([],Plant.Data.Timestamp,Plant.Data.Weather.(S{j}));
            else
                Plant.Data.HistProf.(S{j}) = TypicalDay([],TestData.Timestamp,TestData.Weather.(S{j}));
            end
        end
    end
    if isfield(Plant.Data,'Hydro')
        nodes = length(TestData.Hydro.SourceSink(1,:));
        if isfield(Plant.Data.Hydro,'SourceSink') && nodes == length(Plant.Data.Hydro.SourceSink(1,:))
            for n = 1:1:nodes
                Plant.Data.HistProf.Hydro.SourceSink(n) = {TypicalDay([],Plant.Data.Hydro.Timestamp,Plant.Data.Hydro.SourceSink(:,n))};
            end
        else
            for n = 1:1:nodes
                Plant.Data.HistProf.Hydro.SourceSink(n) = {TypicalDay([],TestData.Hydro.Timestamp,TestData.Hydro.SourceSink(:,n))};
            end
        end
    end
    if strcmp(Plant.optimoptions.forecast,'Surface')
        calculateHistoricalFit; %% calculate surface fits used in forecasting
    end
end

%% Load generators, & build QP matrices
loadGenerator % Loads generators stored in Plant.Generator
loadBuilding %Loads any buildings into QPform
build_subNet
findbuffer %need to wait until everything else has been loaded to load the buffer, because it is reliant on how quickly everything else can charge the system
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
Plant.OpMatA = buildMatrices('A',dt); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B',dt);%build quadratic programming matrices for FitB
Plant.OneStep = buildMatrices1Step;%build quadratic programming matrices for 1Step at interval spacing of dt

if strcmp(Plant.optimoptions.method,'Control')
    A.Horizon = Plant.optimoptions.Resolution;%the horizon is the resolution
    A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = buildTimeVector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        Plant.Online(t) = buildMatrices('OpMatB',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end

