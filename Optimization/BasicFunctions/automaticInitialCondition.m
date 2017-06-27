function IC = automaticInitialCondition(Data_t0)
global Plant DateSim OnOff CurrentState
nG = length(Plant.Generator);
if strcmp(Plant.optimoptions.solver,'NREL')
    %skip initialization
    IC = zeros(1,nG);
else
    scaleCost = updateGeneratorCost(DateSim);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
    Data_t0.Renewable = zeros(1,length(Plant.Generator));
    for i = 1:1:length(Plant.Generator)
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Data_t0.Renewable(i) = RenewableOutput(i,DateSim,'Actual');
        end
    end
    if isfield(Plant.Network,'Hydro') % ~any(Plant.OpMatA.Organize.Dispatchable)
        %need a new initialization
    else
        IC = StepByStepDispatch(Data_t0,scaleCost,Plant.optimoptions.Resolution,[],'',[]);
    end
    for i=1:1:nG
        if isfield(Plant.Generator(i).QPform,'Stor') && Plant.Generator(i).Enabled
            IC(i) = 0.5*Plant.Generator(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
        elseif isfield(Plant.Generator(i).QPform,'Stor')
            IC(i) = 0; %storage that is disabled has 0 IC
        end
    end

    %% specify initial river flow and spillway flows: IC (nL)
    if isfield(Plant.Network,'Hydro')
        NodeNames = cell(length(Plant.Network),1);
        for i = 1:1:length(Plant.Network)
            NodeNames(i) = {Plant.Network(i).name};
        end
        networkNames = fieldnames(Plant.Network);
        networkNames = networkNames(~strcmp('name',networkNames));
        networkNames = networkNames(~strcmp('Equipment',networkNames));
        for net = 1:1:length(networkNames)
            nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
        end
        nLcum = 0; %cumulative line #
        for net = 1:1:length(networkNames)
            if strcmp(networkNames{net},'Hydro') %all hydro lines have initial state
                for i = 1:1:nLinet(net) 
                    %find name of upstream node
                    name = Plant.subNet.lineNames.Hydro{i};
                    r = strfind(name,'_');
                    I = find(strcmp(name(1:r(1)-1),Plant.Data.Hydro.Nodes));
                    if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                        IC(nG+nLcum+i) = Data_t0.Hydro.OutFlow(I);
                    else IC(nG+nLcum+i) = Data_t0.Hydro.SpillFlow(I);
                    end
                end
            end
            nLcum = nLcum+nLinet(net);
        end
    end

    OnOff = true(1,nG);
    for i = 1:1:nG
        if isempty(strfind(Plant.Generator(i).Type,'Storage')) && isempty(strcmp(Plant.Generator(i).Type,'Utility')) && isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
            states = Plant.Generator(i).QPform.states(:,end);
            LB = 0;
            for j = 1:1:length(states);
                LB = LB + Plant.Generator(i).QPform.(states{j}).lb;
            end
            if IC(i)<LB
                OnOff(i) = false;
            end
        end
    end
end

CurrentState.Generators=IC(1:nG);
CurrentState.Lines=IC(nG+1:end);