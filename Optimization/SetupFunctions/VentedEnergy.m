function [eVented,xL,lb,ub] = VentedEnergy(xL,lb,ub,param)
%add states at each node to allow puroseful energy discharge (venting heat or cooling)
global Plant
if strcmp(param,'H') 
    net = 'DistrictHeat';
    opt = 'excessHeat';
elseif strcmp(param,'C')
    net = 'DistrictCool';
    opt = 'excessCool';
end
if any(strcmp(net,fieldnames(Plant.subNet))) && Plant.optimoptions.(opt) == 1
    %%find maximum heat production possible
    maxGen = 0;
    nG = length(Plant.Generator);
    for i = 1:1:nG
        if isfield(Plant.Generator(i).QPform.output,param)
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
                maxGen = maxGen + Plant.Generator(i).QPform.Stor.PeakDisch;
            else
                if strcmp(param,'H') && isfield(Plant.Generator(i).QPform.output,'E')
                    maxGen = maxGen + max(Plant.Generator(i).Size*Plant.Generator(i).Output.Capacity./Plant.Generator(i).Output.Electricity.*Plant.Generator(i).Output.Heat);
                elseif Plant.Generator(i).QPform.output.(param)(end)>0
                    maxGen = maxGen + max(max(Plant.Generator(i).QPform.output.(param)))*Plant.Generator(i).Size;
                end
            end
        end
    end
    %%assume heat can be lost any any node in the network that has a device producing heat
    n = length(Plant.subNet.(net).nodes);
    ventNodes = 0;
    eVented =zeros(n,1); %matrix for the state associated with venting heat at each district heating node, at each time step
    for i = 1:1:n
        genI = Plant.subNet.(net).Equipment{i};%%identify generators at this node
        for j = 1:1:length(genI)
            if isfield(Plant.Generator(genI(j)).QPform.output,'H')
                ventNodes = ventNodes + 1; % Add single state for heat that is ventd to make energy equality true
                eVented(1,i) = (xL+ventNodes);
                break
            end
        end
    end
    xL = xL+ventNodes;
    lb(end+1:end+ventNodes,1) = 0;
    ub(end+1:end+ventNodes,1) = maxGen;
elseif any(strcmp(net,fieldnames(Plant.subNet)))
    eVented = zeros(length(Plant.subNet.(net).nodes),1);
else
    eVented = [];
end
end%Ends function VentedHeat