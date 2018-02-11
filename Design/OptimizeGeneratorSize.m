function OptimizeGeneratorSize
%OPTIMIZEGENERATORSIZE Optimize the size of the selected generator.

global Plant testSystems SYSINDEX GENINDEX

% load plant
Plant = testSystems(SYSINDEX);
Gen0 = testSystems(SYSINDEX).Generator(GENINDEX);
genSize = [(0.2*Gen0.Size : 0.2*Gen0.Size : Gen0.Size),1.5*Gen0.Size,2*Gen0.Size];
handles = guihandles;
if get(handles.DesignDay,'Value') == 1
    DD = true;
else
    DD = false;
end
Years = str2double(get(handles.NPC_Years,'String'));%choose years in GUI
% find net present cost for each size
npc = zeros(size(genSize));
for i = 1:1:length(genSize)
    scale = genSize(i)/Gen0.Size;
    npc(i) = GeneratorNpc(scale,Gen0,DD,Years);
end

% get optimal size
[~,iMinCost] = min(npc);
optimSize = genSize(iMinCost);

% put Plant back into testSystems
scale = optimSize/Gen0.Size;
GenNew = updateComponentSpec(Gen0,'UB',Gen0.Size*scale);
F = fieldnames(GenNew);
for j = 1:1:length(F)
    testSystems(SYSINDEX).Generator(GENINDEX).(F{j}) = GenNew.(F{j});
end
DesignCosts(SYSINDEX,Years,Plant.Costs.Equipment);
testSystems(SYSINDEX).Design = [];
MainScreen1('popupmenu_Callback')
end%Ends function OptimizeGeneratorSize
