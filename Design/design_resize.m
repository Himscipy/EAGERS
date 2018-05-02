function [gen,equip_costs] = design_resize(gen,equip_costs,size,index_size)
% for each generator specified by the input indices, update Plant with its
% newly scaled version
for j = 1:1:length(index_size)
    i = index_size(j);
    scale = size(j) / gen(i).Size;
    gen(i) = update_component_spec(gen(i), 'UB', size(j));   
    equip_costs(i).Cost = equip_costs(i).Cost * scale;
    equip_costs(i).OandM = equip_costs(i).OandM * scale;
end
end%ends function design_resize