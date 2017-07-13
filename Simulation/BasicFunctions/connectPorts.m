function connectPorts
global modelParam Outlet Inlet
CompNames = fieldnames(modelParam.Components);
controls = fieldnames(modelParam.Controls);
list = [CompNames;controls;];
for k = 1:1:length(list)
    block = list{k};
    if any(strcmp(controls,block))
        Co = 'Controls';
    else Co = 'Components';
    end
    list2 = modelParam.(Co).(block).InletPorts;
    for i = 1:1:length(list2)
        port = list2{i};
        if isfield(modelParam.(Co).(block).(port),'IC')
                Inlet.(block).(port) = modelParam.(Co).(block).(port).IC; %use initial condition, update later if connected to an outlet
            else Inlet.(block).(port) = []; %don't have an IC for this inlet port, will get updated later
        end
        if length(modelParam.(Co).(block).connections)<i || isempty(modelParam.(Co).(block).connections{i})
            modelParam.(Co).(block).(port).connected={};
        else
            if ischar(modelParam.(Co).(block).connections{i})
                modelParam.(Co).(block).(port).connected = modelParam.(Co).(block).connections(i);
            else
                modelParam.(Co).(block).(port).IC = modelParam.(Co).(block).connections{i};
                modelParam.(Co).(block).(port).connected={};
            end
        end
    end
    list3 = modelParam.(Co).(block).OutletPorts;
    for i = 1:1:length(list3)
        port = list3{i};
        Outlet.(block).(port) = modelParam.(Co).(block).(port).IC;
    end   
end
end %Ends function connectPorts