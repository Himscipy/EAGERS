function maxUtilitySellback
%identify upper bound for utility states (helps with scaling)
global Plant TestData
nG = length(Plant.Generator);
maxDCgen = 0;
if isfield(Plant.subNet,'DirectCurrent')
    for n = 1:1:length(Plant.subNet.DirectCurrent.nodes)
        equip = Plant.subNet.DirectCurrent.Equipment{n};
        for j = 1:1:length(equip)
            if ~strcmp(Plant.Generator(equip(j)).Type,'AC_DC')
                maxDCgen = maxDCgen+ max(0,Plant.Generator(equip(j)).Size*Plant.Generator(equip(j)).QPform.output.DC(1));
            end
        end
    end
end
maxACgen = 0;
if isfield(Plant.subNet,'Electrical')
    for n = 1:1:length(Plant.subNet.Electrical.nodes)
        equip = Plant.subNet.Electrical.Equipment{n};
        for j = 1:1:length(equip)
            if ~strcmp(Plant.Generator(equip(j)).Type,'AC_DC')
                maxACgen = maxACgen+ max(0,Plant.Generator(equip(j)).Size*Plant.Generator(equip(j)).QPform.output.E(1));
            end
        end
    end
end
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Utility') && ~isempty(Plant.Generator(i).QPform.states)%%avoid things like gas utility with no states
        %identify the network
        out = char(fieldnames(Plant.Generator(i).QPform.output));
        if strcmp(out,'E')
            net = 'Electrical';
        elseif strcmp(out,'H')
            net = 'DistrictHeat';
        elseif strcmp(out,'C')
            net = 'DistrictCool';
        elseif strcmp(out,'W')
            net = 'Hydro';
        end
        if isfield(TestData,'Demand')
            Plant.Generator(i).QPform.X.ub = 10*max(TestData.Demand.(out));%max Purchase
        else
            Plant.Generator(i).QPform.X.ub = 1e6; %arbitrary upper bound that is not inf
        end
        if strcmp(net,'Electrical') && isfield(Plant.Generator(i).QPform,'Y')
            Plant.Generator(i).QPform.Y.ub = maxACgen + maxDCgen;
        end
    end
    if strcmp(Plant.Generator(i).Type,'AC_DC')
        Plant.Generator(i).QPform.A.ub = maxACgen;
        Plant.Generator(i).QPform.B.ub = maxDCgen;
    end
end 