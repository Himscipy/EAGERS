global Plant RestartTime GenAvailTime DateSim
nG = length(Plant.Generator);
if any(strcmp(Plant.optimoptions.method,{'Control'}))
    for i = 1:1:nG
        if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
            RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
        else
            RestartTime(i) = 0;
        end
    end
    GenAvailTime = ones(1,nG).*DateSim; %  Global variable needed in controller mode
end