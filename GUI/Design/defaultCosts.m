function Costs = defaultCosts(Costs,Gen)
%DEFAULTCOSTS
if isempty(Costs) || ~isfield(Costs,'DiscountRate')
    Costs.DiscountRate = 2;% implies a 2% discount rate for NPC calculations
end
nG = length(Gen);
for i = 1:1:nG
    Type = Gen(i).Type;
    if ~strcmp(Type,'Utility') && (~isfield(Costs,'Equipment') || length(Costs.Equipment)<i || isempty(Costs.Equipment(i)))
        Costs.Equipment(i).Name = Gen(i).Name;
        if (strcmp(Type,'CHP Generator') || strcmp(Type,'Electric Generator')) && Gen(i).VariableStruct.isFuelCell
            costPerkW = 3000;
            OM = 100;
        elseif strcmp(Type,'CHP Generator') || strcmp(Type,'Electric Generator')
            costPerkW = 1000;
            OM = 50;
        elseif strcmp(Type,'Chiller') %chillers
            if strcmp(Gen(i).Source,'Electricity')
                % electric chiller
                costPerkW = 100;
                OM = 10;
            else
                % absorption chiller
                costPerkW = 200;
                OM = 20;
            end
        elseif strcmp(Type,'Heater') 
            costPerkW = 100;
            OM = 3;
        elseif strcmp(Type,'Thermal Storage') 
            costPerkW = 20;
            OM = 1;
        elseif strcmp(Type,'Electric Storage')
            costPerkW = 500;
            OM = 10;
        elseif strcmp(Type,'Solar')
            costPerkW = 500;
            OM = 5;
        end
        Costs.Equipment(i).Cost = costPerkW*Gen(i).Size;
        Costs.Equipment(i).OandM = OM*Gen(i).Size;
        Costs.Equipment(i).Financed = 100;
        Costs.Equipment(i).LoanRate = 6;
        Costs.Equipment(i).LoanTerm = 15;
    end
end