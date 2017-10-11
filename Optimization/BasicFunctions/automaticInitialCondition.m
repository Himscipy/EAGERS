function IC = automaticInitialCondition(Data_t0)
global Plant DateSim OnOff CurrentState
if strcmp(Plant.optimoptions.solver,'NREL')
    nG = length(Plant.Generator);%skip initialization
    IC = zeros(1,nG);
    
    for i=1:1:nG
        if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
            IC(i) = 0.5*Plant.Generator(i).Size*(Plant.Generator(i).VariableStruct.MaxDOD/100); % IC = halfway charged energy storage
        end
    end
    CurrentState.Generators=IC(1:nG);
else
    nG = length(Plant.Generator);
    nB = length(Plant.Building);
    nL = length(Plant.OpMatA.Organize.IC)-nG-nB;
    scaleCost = updateGeneratorCost(DateSim);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
    Data_t0.Renewable = zeros(1,length(Plant.Generator));
    for i = 1:1:nG
        if ismember(Plant.Generator(i).Type,{'Solar';'Wind';}) && Plant.Generator(i).Enabled
            Data_t0.Renewable(i) = RenewableOutput(i,DateSim,'Actual');
        end
    end
    %% Specify initial building temperatures
    CurrentState.Buildings = [];
    for i = 1:nB
        CurrentState.Buildings(i) = 21.1; %initial temperature guess
    end
    
    IC = StepByStepDispatch(Data_t0,scaleCost,Plant.optimoptions.Resolution,'',[]);
    OnOff = true(1,nG);
    for i=1:1:nG
        if ismember(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';})
            IC(i) = 0.5*Plant.Generator(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
        elseif Plant.OpMatA.Organize.Dispatchable(i)
            states = Plant.Generator(i).QPform.states(:,end);
            LB = 0;
            for j = 1:1:length(states)
                LB = LB + Plant.Generator(i).QPform.(states{j}).lb;
            end
            if IC(i)<LB
                OnOff(i) = false;
            end
        end
    end
    if isfield(Plant.Network,'Hydro')
        CurrentState.Hydro = zeros(1,nG);%SOC for reserviors,  IC is the initial power production
        for n=1:1:length(Plant.subNet.Hydro.nodes)
            equip = Plant.subNet.Hydro.Equipment{n};
            for k = 1:1:length(equip)
                i = equip(k);
                if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                    CurrentState.Hydro(i) = 0.5*Plant.Generator(i).QPform.Stor.UsableSize; % IC = halfway charged energy storage
%                     IC(i) = (Data_t0.Hydro.OutFlow(n) - Data_t0.Hydro.SpillFlow(n))/Plant.Generator(i).QPform.Stor.Power2Flow; %inital power
                end
            end
        end
        %load the previous 24 hours of SourceSink & Outflows
        interval = Plant.Data.Hydro.Timestamp(2)-Plant.Data.Hydro.Timestamp(1);
        if Plant.Data.Hydro.Timestamp(1)>(DateSim-1)
            %take 1st 24 hours and treat as yesterday
            D = DateSim;
            while Plant.Data.Hydro.Timestamp(1)>D
                D = D+1;
            end
            xi = nnz(Plant.Data.Hydro.Timestamp<=D);
            xf = nnz(Plant.Data.Hydro.Timestamp<D+1+interval);
        else %take last 24 hours
            xi = nnz(Plant.Data.Hydro.Timestamp<=(DateSim-1));
            xf = nnz(Plant.Data.Hydro.Timestamp<DateSim+interval);
        end
        n = length(Plant.subNet.Hydro.nodes);
        Plant.Data.HydroHistory.Timestamp = linspace((DateSim-1),DateSim,1/interval+1)';
        Plant.Data.HydroHistory.SourceSink = Plant.Data.Hydro.SourceSink(xi:xf,1:n);
        Plant.Data.HydroHistory.OutFlow = Plant.Data.Hydro.OutFlow(xi:xf,1:n);
    end
%     %% specify initial river flow : IC (nL)
%     networkNames = fieldnames(Plant.subNet);
%     if ismember('Hydro',networkNames) %all hydro lines have initial state
%         for i = 1:1:length(Plant.subNet.Hydro.lineNames)
%             %find name of upstream node
%             name = Plant.subNet.Hydro.lineNames{i};
%             r = strfind(name,'_');
%             node = nonzeros((1:1:length(Plant.Data.Hydro.Nodes))'.*strcmp(name(1:r(1)-1),Plant.Data.Hydro.Nodes));
%             IC(nG+Plant.subNet.Hydro.lineNumber(i)) = Data_t0.Hydro.OutFlow(node);
%         end
%     end

    CurrentState.Generators=IC(1:nG);%Moved this from line 43
    CurrentState.Lines= [];%IC(nG+1:nG+nL); %do we need this anymore?
    IC(nG+nL+1:nG+nL+nB) = CurrentState.Buildings;
end

end%Ends function automaticInitialCondition