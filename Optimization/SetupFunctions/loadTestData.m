global Plant TestData
if isempty(TestData) 
    if~isempty(Plant.Data)
        TestData = Plant.Data;%% Revise this so you can pull from more than what is loaded in Plant
        if isfield(TestData,'HistProf')
            TestData = rmfield(TestData,'HistProf');
        end
        if isfield(Plant.Data,'Hydro')
            TestData = rmfield(TestData,'Hydro');
            TestData.Hydro.SourceSink = Plant.Data.Hydro.SourceSink;
            TestData.Hydro.OutFlow = Plant.Data.Hydro.OutFlow;
        end            
    end
    if ~isfield(TestData,'Timestamp')
        nS = 365*24/Plant.optimoptions.Resolution+1;
        TestData.Timestamp = linspace(datenum([2017,1,1]),datenum([2018,1,1]),nS)';
    end
    S = {'Tdb';'Twb';'irradDireNorm';};
    for j = 1:1:length(S)
        if ~isfield(TestData,'Weather') || ~isfield(TestData.Weather,S{j})
            D = datevec(TestData.Timestamp(1));
            TestData.Weather.(S{j}) = interp1(linspace(0,8760,length(Plant.Weather.Tdb)+1)',[Plant.Weather.(S{j})(1); Plant.Weather.(S{j})],mod(24*(TestData.Timestamp - datenum([D(1),1,1])),8760)); 
        end
    end

    if ~isempty(Plant.Building)
        %% add building data to extra columns in test data, or do nothing here?
        if ~isfield(TestData,'Demand')
            nS = length(TestData.Timestamp);
            TestData.Demand.E = zeros(nS,0);
            TestData.Demand.C = zeros(nS,0);
            TestData.Demand.H = zeros(nS,0);
        end
        for i = 1:1:length(Plant.Building)
            %use the network structure to find the location of the building
            Location = [];
            [Equipment,InteriorLighting,ExteriorLighting,Cooling,Heating,FanPower,OtherLoads] = BuildingProfile(Plant.Building(i),Plant.Weather,TestData.Timestamp,Location);
            % Compare2Eplus(Plant.Building(i),Plant.Weather,TestData.Timestamp);
            %% need conditions of if there are chillers and heaters only for electric
            Electric = Equipment + InteriorLighting + ExteriorLighting + FanPower + Cooling/Plant.Building(i).VariableStruct.COP_C + Heating/Plant.Building(i).VariableStruct.COP_H + OtherLoads;
            TestData.Demand.E(:,end+1) = Electric;
            TestData.Demand.C(:,end+1) = Cooling;
            TestData.Demand.H(:,end+1) = Heating;
        end
    end
end