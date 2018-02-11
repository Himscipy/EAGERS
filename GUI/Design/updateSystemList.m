global testSystems SYSINDEX
nG = length(testSystems(SYSINDEX).Generator);
list = cell(nG,1);
for i = 1:1:nG
    list(i) = {strcat(testSystems(SYSINDEX).Generator(i).Type,'.',testSystems(SYSINDEX).Generator(i).Name)};
end
testSystems(SYSINDEX).Network.Equipment = list;
