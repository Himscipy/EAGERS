function NumSteps = preAllocateSpace(mode)
global Plant TestData
nG = length(Plant.Generator);
if isfield(Plant,'Building') && ~isempty(Plant.Building)
    nB = length(Plant.Building);
else
    nB = 0;
end
[~,n] = size(Plant.OneStep.organize);
nL = n-nG-nB;
if isfield(TestData,'Hydro')
    nD = length(TestData.Hydro.SourceSink(1,:));
else
    nD = 0;
end
nS = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
F = fieldnames(TestData);
F = F(~strcmp('Timestamp',F));
F = F(~strcmp('RealTimeData',F));
if strcmp(mode,'Design')
    Plant.Design.Timestamp = zeros(NumSteps,1);   
    Plant.Design.GeneratorState = zeros(NumSteps,nG);
    Plant.Design.LineFlows = zeros(NumSteps,nL);
    Plant.Design.Buildings = zeros(NumSteps,nB);
    Plant.Design.hydroSOC = zeros(NumSteps,nD);
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                if isnumeric(TestData.(F{j}).(S{i}))
                    l = length(TestData.(F{j}).(S{i})(1,:));
                    Plant.Design.(F{j}).(S{i}) = zeros(NumSteps,l);
                end
            end
        elseif ~isempty(TestData.(F{j})) && isnumeric(TestData.(F{j}))
            l = length(TestData.(F{j})(1,:));
            Plant.Design.(F{j}) = zeros(NumSteps,l);
        end
    end
elseif strcmp(mode,'Dispatch')
    Plant.Dispatch.Timestamp = zeros(NumSteps,1);
    Plant.Dispatch.GeneratorState = zeros(NumSteps,nG);
    Plant.Dispatch.LineFlows = zeros(NumSteps,nL);
    Plant.Dispatch.Buildings = zeros(NumSteps,nB);
    Plant.Dispatch.hydroSOC = zeros(NumSteps,nD);
    Plant.Predicted.GenDisp = zeros(nS,nG,NumSteps);
    Plant.Predicted.LineFlows = zeros(nS,nL,NumSteps);
    Plant.Predicted.Buildings = zeros(nS,nB,NumSteps);
    Plant.Predicted.hydroSOC = zeros(nS,nD,NumSteps);
    Plant.Predicted.Timestamp = zeros(nS,NumSteps);
    Plant.Predicted.Cost = zeros(NumSteps,1);
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                if isnumeric(TestData.(F{j}).(S{i}))
                    l = length(TestData.(F{j}).(S{i})(1,:));
                    Plant.Predicted.(F{j}).(S{i}) = zeros(NumSteps,nS,l);
                    Plant.Dispatch.(F{j}).(S{i}) = zeros(NumSteps,l);
                end
            end
        elseif ~isempty(TestData.(F{j})) && isnumeric(TestData.(F{j}))
            l = length(TestData.(F{j})(1,:));
            Plant.Predicted.(F{j}) = zeros(NumSteps,nS,l);
            Plant.Dispatch.(F{j}) = zeros(NumSteps,l);
        end
    end
    Plant.RunData = Plant.Dispatch;
end
end%Ends preAllocateSpace
