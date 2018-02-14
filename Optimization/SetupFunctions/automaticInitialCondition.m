function automaticInitialCondition(Data_t0)
global Plant DateSim OnOff CurrentState TestData
CurrentState=[];
nG = length(Plant.Generator);
CurrentState.Generators = zeros(1,nG);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
    Data_t0.Building.ExternalGains = ForecastExternalGains(Data_t0);
    Data_t0.Building = ForecastBuilding(Data_t0.Weather,Data_t0.Timestamp,Data_t0.Building);
else
    nB = 0;
    CurrentState.Buildings = zeros(2,nB);
end
nL = length(Plant.OpMatA.Organize.IC)-nG-nB;
CurrentState.Lines = zeros(1,nL);
scaleCost = updateGeneratorCost(DateSim,Plant.Generator);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
Data_t0.Renewable = zeros(1,length(Plant.Generator));
for i = 1:1:nG
    if ismember(Plant.Generator(i).Type,{'Solar';'Wind';}) && Plant.Generator(i).Enabled
        Data_t0.Renewable = RenewableOutput(DateSim,Data_t0.Weather.irradDireNorm);
    end
end
IC = StepByStepDispatch(Data_t0,scaleCost,Plant.optimoptions.Resolution,'',[]);

LB = zeros(1,nG);
for i=1:1:nG
    if Plant.OpMatA.Organize.Dispatchable(i)
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            LB(i) = LB(i) + Plant.Generator(i).QPform.(states{j}).lb(2);
        end
    elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
        IC(i) = 0.5*Plant.Generator(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
OnOff = IC>LB;

if isfield(Plant.Network,'Hydro')
    CurrentState.Hydro = zeros(1,nG);%SOC for reserviors,  IC is the initial power production
    for n=1:1:length(Plant.subNet.Hydro.nodes)
        equip = Plant.subNet.Hydro.Equipment{n};
        for k = 1:1:length(equip)
            i = equip(k);
            if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                j = datevec(TestData.Hydro.Timestamp(1));
                j = j(2);
                CurrentState.Hydro(i) = Plant.Generator(i).VariableStruct.StartingPoint(j)*Plant.Generator(i).QPform.Stor.UsableSize;
            end 
        end 
    end 
end
CurrentState.Generators=IC(1:nG);
CurrentState.Lines= IC(nG+1:nG+nL); %do we need this anymore?
end%Ends function automaticInitialCondition