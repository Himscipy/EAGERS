%% calculate the upper and lower buffer thresholds for each storage system
global Plant
nG = length(Plant.Generator);
BuffPerc = Plant.optimoptions.Buffer; % percentage for buffer on storage
genNames = cell(nG,1);
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
for net = 1:1:length(networkNames)
    %identify equipment at each subNet node
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
        include = {'Electric Generator','CHP Generator'};
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
        include = {'Heater';'CHP Generator';};
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
        include = {'Chiller'}; %switched this with DistrictCooling include
    elseif strcmp(networkNames{net},'Hydro')
        out = 'W';
        include = {}; % included this to get past line 31
    end
    for m = 1:1:length(Plant.subNet.(networkNames{net}))
        equip = Plant.subNet.(networkNames{net})(m).Equipment;
        chargeCapacity = 0; 
        for j = 1:1:length(equip)
            if ismember(Plant.Generator(equip(j)).Type,include)
                if strcmp(out,'H') && strcmp(Plant.Generator(equip(j)).Type,'CHP Generator')
                    chargeCapacity = chargeCapacity+Plant.Generator(equip(j)).Size*Plant.Generator(equip(j)).QPform.output.H/Plant.Generator(equip(j)).QPform.output.E; %CHP generators heat ratio
                else chargeCapacity = chargeCapacity+Plant.Generator(equip(j)).Size;
                end
            end
        end
        for j = 1:1:length(equip)
            if isfield(Plant.Generator(equip(j)).QPform,'Stor') && any(strcmp(Plant.Generator(equip(j)).QPform.states,'W'))%storage with a buffer state
                if any(strcmp('W',out))
                    %hydro
                    dischargeCapacity = (Plant.Generator(equip(j)).VariableStruct.MaxGenFlow + Plant.Generator(equip(j)).VariableStruct.MaxSpillFlow)/12.1; %flow rate in 1000 ft^3 converted to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                    dischargeCapacity = dischargeCapacity*Plant.optimoptions.Horizon/10; %amount the resevoir can discharge in 10% of the dispatch horizon.
                    Buffer = min((BuffPerc/100)*Plant.Generator(equip(j)).QPform.Stor.UsableSize,dischargeCapacity);
                else
                    chargeCapacity = chargeCapacity*2*Plant.optimoptions.Resolution; %amount the local systems can charge the storage in two steps of the dispatch horizon.
                    Buffer = min((BuffPerc/100)*Plant.Generator(equip(j)).QPform.Stor.UsableSize,chargeCapacity);
                end
                Plant.Generator(equip(j)).QPform.link.bineq(end-1) = -Buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
                Plant.Generator(equip(j)).QPform.link.bineq(end) = Plant.Generator(equip(j)).QPform.Stor.UsableSize-Buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
                Plant.Generator(equip(j)).QPform.Z.ub = Buffer;
                Plant.Generator(equip(j)).QPform.W.ub = Buffer;
            end
        end
    end
end