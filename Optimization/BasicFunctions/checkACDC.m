function nG = checkACDC(mode)
global testSystems SYSINDEX Model_dir TestData Plant
%double check & add AC_DC conversion if neccessary
hasDC = false;
hasAC = false;
hasAC_DC = false;
needAC = false;
needDC = false;
if strcmp(mode,'plan')
    GenList = testSystems(SYSINDEX).Generator;
elseif strcmp(mode,'control')
    GenList = Plant.Generator;
end
nG = length(GenList);
for i = 1:1:nG
    if strcmp(GenList(i).Type,'AC_DC') 
        hasAC_DC = true;
    end
    if isfield(GenList(i).Output,'DirectCurrent')
        hasDC = true;
    end
    if isfield(GenList(i).Output,'Electricity') || strcmp(GenList(i).Type,'Hydro Storage')
        hasAC = true;
    end
end
if isfield(TestData,'Demand')
    if isfield(TestData.Demand,'E')
        needAC = true;
    end
    if isfield(TestData.Demand,'DC')
        needDC = true;
    end
else
    needAC = true;
    if isfield(TestData.Building,'DCloads')
        needDC = true;
    end
end
if ~hasAC_DC && (needAC && ~hasAC) || (needDC && ~hasDC)
    load(fullfile(Model_dir,'System Library','AC_DC','IdealConverter.mat'))
    GenList(nG+1) = component;
    nG = nG+1;
    if strcmp(mode,'plan')
        testSystems(SYSINDEX).Generator = GenList;
    elseif strcmp(mode,'control')
        Plant.Generator = GenList;
    end
end
end%ends function checkACDC