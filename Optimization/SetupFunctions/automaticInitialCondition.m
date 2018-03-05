function automaticInitialCondition(Data_t0)
global Plant DateSim OnOff CurrentState
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    Data_t0.Building.ExternalGains = zeros(1,length(Plant.Building));
    for i = 1:1:length(Plant.Building)
        Location = Plant.subNet.Electrical.Location(Plant.Building(i).QPform.Electrical.subnetNode);
        SG = SolarGain(Plant.Building(i),DateSim,Location,Data_t0.Weather);
        Data_t0.Building.ExternalGains(1,i) = SG.Walls + SG.Roof;
    end
    Data_t0.Building = ForecastBuilding(Data_t0.Weather,Data_t0.Timestamp,Data_t0.Building);
end
Type = cell(nG,1);
scaleCost = updateGeneratorCost(DateSim,Plant.Generator);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
for i = 1:1:nG
    Type(i) = {Plant.Generator(i).Type};
end
if any(strcmp(Type,'Solar')) || any(strcmp(Type,'Wind'))
    Data_t0.Renewable = RenewableOutput(DateSim,Data_t0.Weather.irradDireNorm);
end
IC = StepByStepDispatch(Data_t0,scaleCost,Plant.optimoptions.Resolution,[]);

LB = zeros(1,nG);
for i=1:1:nG
    if Plant.OpMatA.Organize.Dispatchable(i)
        states = Plant.Generator(i).QPform.states(1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end))),end);
        for j = 1:1:length(states)
            LB(i) = LB(i) + Plant.Generator(i).QPform.(states{j}).lb(2);
        end
    elseif ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';'Hydrogen Storage';})
        IC(i) = 0.5*Plant.Generator(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
OnOff = IC(1:nG)>LB;
CurrentState.Generators=IC(1:nG);
if isfield(Plant.Network,'Hydro') && isfield(Plant,'WYForecast')
    nY = length(Plant.WYForecast);
    prevForcast = Plant.WYForecast{nY};
    for i = 1:1:length(prevForcast.hydroSOC(1,:))
        CurrentState.Hydro(i) = interp1(prevForcast.Timestamp,prevForcast.hydroSOC(:,i),DateSim);
    end
end
end%Ends function automaticInitialCondition