function hydro_diffTime
global Plant
%diffTime = upriver(dam) time from columbia rmouth - downriver(dam) time from columbiar mouth

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));

nH = length(Plant.Generator);
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Hydro')
        for dam = 1:1:nH
            Plant.Generator(dam).VariableStruct.diffTime = [];
            for i = 1:1:nH
                if nnz(find(strcmp(strcat(Plant.Generator(i).Name,'_',networkNames{net},'_',Plant.Generator(dam).Name),Plant.subNet.lineNames.(networkNames{net})),1,'first'));
                    Plant.Generator(dam).VariableStruct.diffTime(end+1)= Plant.Generator(i).VariableStruct.transTime - Plant.Generator(dam).VariableStruct.transTime;
                end 
            end  
        end 
    end 
end
