function manualInitialCondition
global Plant OnOff CurrentState
CurrentState = [];
nG = length(Plant.Generator);
nB = length(Plant.Building);
nL = length(Plant.OpMatA.Organize.IC)-nG-nB;
list = {};
list2 = {};
sizes = {};
Index =[];
IC = zeros(1,nG+nL);
CurrentState.Generators = zeros(1,nG);
CurrentState.Buildings = zeros(2,nB)+17.4;
include = {'CHP Generator', 'Electric Generator','Chiller','Heater','Hydro Storage'};
for i = 1:1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        list(end+1) = {strcat(Plant.Generator(i).Name,' --- Electric Output')};
        list2(end+1) = {strcat(Plant.Generator(i).Name,' --- % of max capacity')};
        sizes(end+1) = {'50'};
        Index(end+1) = i;
    elseif isfield(Plant.Generator(i).OpMatA,'Stor')
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
list = [list,list2];
sizes = [sizes,50*ones(1,length(list2))];
Input = str2double(inputdlg(list,'Specify Initial Condition (kW) or State of Charge (%)',1,sizes));
IC(Index) = Input(1:length(Index));
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor') && ~strcmp(Plant.Generator(i).Type,'Hydro Storage')
        IC(i) = IC(i)/100*Plant.Generator(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
%% specify initial river flow and spillway flows: IC (nL)
networkNames = fieldnames(Plant.subNet);
if any(strcmp('Hydro',networkNames))
    for n = 1:1:length(Plant.subNet.Hydro.nodes) 
        i = Plant.subNet.Hydro.Equipment{n};
        if strcmp(Plant.Generator(i).Type,'Hydro Storage')
            CurrentState.Hydro(n) =  Input(length(Index)+n);
        end
    end
end
OnOff = true(1,nG);
for i = 1:1:nG
    if isempty(strfind(Plant.Generator(i).Type,'Storage')) && isempty(strfind(Plant.Generator(i).Type,'Utility'))
        states = Plant.Generator(i).OpMatB.states;
        LB = 0;
        for j = 1:1:length(states)
            LB = LB + Plant.Generator(i).OpMatB.(states{j}).lb;
        end
        if IC(i)<LB
            OnOff(i) = false;
        end
    end
end
CurrentState.Generators=IC(1:nG);
CurrentState.Buildings(:,1:nB) = 20; %initial temperature guess
end%Ends function manualInitialCondition