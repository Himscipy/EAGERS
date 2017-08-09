function IC = manualInitialCondition
global Plant DateSim OnOff CurrentState
nG = length(Plant.Generator);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG;
list = {};
sizes = {};
Index =[];
IC = zeros(1,nG+nL);
include = {'CHP Generator', 'Electric Generator', 'Chiller','Heater'};
for i = 1:1:length(Plant.Generator)
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        list(end+1) = {strcat(Plant.Generator(i).Name,' --- % of max capacity')};
        sizes(end+1) = {'50'};
        Index(end+1) = i;
    elseif ismember(cellstr(Plant.Generator(i).Type),include)
        list(end+1) = {strcat(Plant.Generator(i).Name,' --- If greater than 0, IC must be between lower bound(',num2str(LB(i)),') and upper bound(',num2str(Plant.Generator(i).Size),').')};
        sizes(end+1) = {'0'};
        Index(end+1) = i;
    elseif ~isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
        %renewable
    elseif ~isempty(strcmp(Plant.Generator(i).Type,'Utility'))
        %utility
     end
end
IC(Index) = str2double(inputdlg(list,'Specify Initial Condition (kW) or State of Charge (%)',1,sizes));
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        IC(i) = IC(i)/100*Plant.Generator(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
%% specify initial river flow and spillway flows: IC (nL)
networkNames = fieldnames(Plant.subNet);
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end

if any(strcmp('Hydro',networkNames))
    for n = 1:1:length(Plant.subNet.Hydro.nodes) 
        name = Plant.subNet.Hydro.lineNames{n};
        r = strfind(name,'_');
        name = name(1:r(1)-1);
        IC(nG+Plant.subNet.Hydro.lineNumber(n)) = str2double(inputdlg(name,'Specify Initial Flow (1000 ft^3/s)',1,0));
        equip = Plant.subNet.Hydro.Equipment{n};
        for k = 1:1:length(equip)
            i = equip(k);
            if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                CurrentState.Hydro(i) =  str2double(inputdlg(Plant.Generator(i).Name,'Specify Initial State of Charge (%)',1,0));
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
OnOff = true(1,nG);
for i = 1:1:nG
    if isempty(strfind(Plant.Generator(i).Type,'Storage')) && isempty(strfind(Plant.Generator(i).Type,'Utility'))
        states = Plant.Generator(i).OpMatB.states;
        LB = 0;
        for j = 1:1:length(states);
            LB = LB + Plant.Generator(i).OpMatB.(states{j}).lb;
        end
        if IC(i)<LB
            OnOff(i) = false;
        end
    end
end

CurrentState.Generators=IC(1:nG);
CurrentState.Lines=IC(nG+1:end);
end%Ends function manualInitialCondition