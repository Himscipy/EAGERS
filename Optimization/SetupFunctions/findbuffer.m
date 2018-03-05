%% calculate the upper and lower buffer thresholds for each storage system
global Plant
nG = length(Plant.Generator);
genNames = cell(nG,1);
networkNames = fieldnames(Plant.subNet);
for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
for net = 1:1:length(networkNames)
    %identify equipment at each subNet node
    for m = 1:1:length(Plant.subNet.(networkNames{net}).nodes)
        equip = Plant.subNet.(networkNames{net}).Equipment{m};
        for j = 1:1:length(equip)
            if isfield(Plant.Generator(equip(j)).QPform,'Stor') && any(strcmp(Plant.Generator(equip(j)).QPform.states,'U'))%storage with a buffer state
                if isfield(Plant.Generator(equip(j)).VariableStruct,'Buffer')
                    BuffPerc = Plant.Generator(equip(j)).VariableStruct.Buffer;% percentage for buffer on storage
                else
                    Plant.Generator(equip(j)).VariableStruct.Buffer = 0;
                    BuffPerc = 0;
                end
                if strcmp(networkNames{net},'Hydro')
                    %hydro
                    dischargeCapacity = (Plant.Generator(equip(j)).VariableStruct.MaxGenFlow + Plant.Generator(equip(j)).VariableStruct.MaxSpillFlow)/12.1; %flow rate in 1000 ft^3 converted to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                    dischargeCapacity = dischargeCapacity*Plant.optimoptions.Horizon/10; %amount the resevoir can discharge in 10% of the dispatch horizon.
                    Buffer = min((BuffPerc/100)*Plant.Generator(equip(j)).QPform.Stor.UsableSize,dischargeCapacity);
                else
                    Buffer = (BuffPerc/100)*Plant.Generator(equip(j)).QPform.Stor.UsableSize;
                end
                Plant.Generator(equip(j)).QPform.link.bineq(end-1) = -Buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
                Plant.Generator(equip(j)).QPform.link.bineq(end) = Plant.Generator(equip(j)).QPform.Stor.UsableSize-Buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
                Plant.Generator(equip(j)).QPform.U.ub = Buffer;
                Plant.Generator(equip(j)).QPform.L.ub = Buffer;
            elseif isfield(Plant.Generator(equip(j)).QPform,'Stor')
                Plant.Generator(equip(j)).VariableStruct.Buffer = 0;
            end
        end
    end
end