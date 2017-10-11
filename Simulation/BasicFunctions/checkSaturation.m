function Inlet = checkSaturation(Inlet,block)
%Check to avoid exceeding bounds on inlets
%This avoids valves in negative position or negative flows
ports = block.InletPorts;
for i = 1:1:length(ports)
    if isfield(block.(ports{i}),'Saturation')
        Sat = block.(ports{i}).Saturation;
    else Sat = [0,inf];
    end
    if isstruct(Inlet.(ports{i}))
        F = fieldnames(Inlet.(ports{i}));
        for j = 1:1:length(F)
            if ~isinf(Sat(1))
                Inlet.(ports{i}).(F{j}) = max(Inlet.(ports{i}).(F{j}),Sat(1));
            end
            if ~isinf(Sat(2))
                Inlet.(ports{i}).(F{j}) = min(Inlet.(ports{i}).(F{j}),Sat(2));
            end
        end
    else
        if ~isinf(Sat(1))
            Inlet.(ports{i}) = max(Inlet.(ports{i}),Sat(1));
        end
        if ~isinf(Sat(2))
            Inlet.(ports{i}) = min(Inlet.(ports{i}),Sat(2));
        end
    end
end
end%ends function checkSaturation