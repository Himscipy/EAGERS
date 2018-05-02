function temperature = cooling_tower_simulate(gen,cool_tower,equip,dispatch)
%need to update with a model of the cooling tower water loop that is not
%precisely the same as in the optimization (break into smaller time steps)
imbalance = 0;
for k = 1:1:length(equip)
    j = equip(k);
    if strcmp(gen(j).Type,'Chiller')
        imbalance = imbalance + dispatch(j) + chill_input(gen(j).QPform,dispatch(j));
    elseif strcmp(gen(j).Type,'Cooling Tower')
        imbalance = imbalance - best_dispatch(j);
    end
end
capacitance = cool_tower.fluid_capacity*cool_tower.fluid_capacitance; %Water capacity in kg and thermal capacitance in kJ/kg*K to get kJ/K
temperature = cool_tower.fluid_temperature + dt*3600/capacitance*imbalance;
end%ends function cooling_tower_simulate

function input = chill_input(qp_form,output)
net_out = 0;
input = 0;
i = 1;
while net_out<output
    seg = min(qp_form.(qp_form.states{i}).ub,output-net_out);
    net_out = net_out + seg;
    if isfield(qp_form.output,'H')
        input = input + seg*qp_form.output.H(i,2);
    end
    if isfield(qp_form.output,'E')
        if length(qp_form.output.E(1,:))>1
            input = input + seg*qp_form.output.E(i,2);
        else
            input = input + seg*qp_form.output.E;
        end
    end
    i = i+1;
end    
end%ends function gen_input