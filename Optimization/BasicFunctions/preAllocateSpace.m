function preAllocateSpace(mode)
global Plant RealTimeData
nG = length(Plant.Generator);
nB = length(Plant.Building);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG-nB;
nS = round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution); %number of steps per dispatch
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1; %number of simulation steps
dems = fieldnames(Plant.Data.Demand);
if strcmp(mode,'Design')
    Plant.Design.Timestamp = zeros(NumSteps,1);
    weath = {'Tdb';'Twb';'irradDireNorm';};
    Plant.Predicted.Weather = [];
    for i = 1:1:length(weath)
        Plant.Design.Weather.(weath{i}) = zeros(NumSteps,1);
    end    
    for i = 1:1:length(dems)
        loads = length(RealTimeData.Demand.(dems{i})(1,:));
        Plant.Design.Demand.(dems{i}) = zeros(NumSteps,loads);
    end
    Plant.Design.GeneratorState = zeros(NumSteps,nG+nL+nB);
elseif strcmp(mode,'Dispatch')
    Plant.Dispatch.Timestamp = zeros(NumSteps,1);
    Plant.Dispatch.GeneratorState = zeros(NumSteps,nG+nL+nB);
    Plant.Predicted.GenDisp = zeros(nS,nG+nL+nB,NumSteps);
    Plant.Predicted.Timestamp = zeros(NumSteps,nS);
    Plant.Predicted.Cost = zeros(NumSteps,1);
    Plant.Predicted.Demand = [];
    for i = 1:1:length(dems)
        loads = length(RealTimeData.Demand.(dems{i})(1,:));
        Plant.Predicted.Demand.(dems{i}) = zeros(NumSteps,nS,loads);
        Plant.Dispatch.Demand.(dems{i}) = zeros(NumSteps,loads);
        Plant.RunData.Demand.(dems{i}) = zeros(NumSteps,loads);
    end
    weath = {'Tdb';'Twb';'irradDireNorm';};
    Plant.Predicted.Weather = [];
    for i = 1:1:length(weath)
        Plant.Predicted.Weather.(weath{i}) = zeros(NumSteps,nS);
    end
    for i = 1:1:nB
        Plant.Dispatch.Building(i).Temperature = zeros(NumSteps,1);
        Plant.RunData.Building(i).Temperature = zeros(NumSteps,1);
    end
end
end%Ends preAllocateSpace
