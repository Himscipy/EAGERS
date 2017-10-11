function [GenDisp,Time,tsim] = DispatchLoop(Date,IC,Forecast,Renewable,PredictDispatch)
%% this dispatch loop trains and then uses an ANN
global Plant

Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
nG = length(Plant.Generator);
nS = length(Time);
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
if isempty(PredictDispatch)
    PredictDispatch = ones(length(Time)+1,1)*IC;
end
scaleCost = updateGeneratorCost(Time/24 + Date); %% All feedstock costs were assumed to be 1 when building matrices 
marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
QP = updateMatrices(Plant.OpMatA,IC,Time,scaleCost,marginCost,Forecast,Renewable,[]);
%% Step 1 Determine initial dispatch
Locked = true(nS+1,nG);
gens = [];
renews = [];
for i = 1:1:nG
    if ~Plant.Generator(i).Enabled
        Locked(:,i) = 0;
    end
    if strcmp(Plant.Generator(i).Type,'Chiller') || strcmp(Plant.Generator(i).Type,'CHP Generator') || strcmp(Plant.Generator(i).Type,'Electric Generator')
        gens(end+1) = i;
    elseif strcmp(Plant.Generator(i).Source, 'Renewable')
        renews(end+1) = i;
    end
end
tic
[FirstDisp,~,Feasible] = DispatchQP(QP,Locked);
marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
renewgen = zeros(nS,nG);
for i = 1:1:length(renews)
    renewgen(:,renews(i)) = RenewableOutput(Plant.Generator(renews(i)).VariableStruct,Date,Time,'Predict');
end
tsim(1,1) = toc;
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    FirstDisp= PredictDispatch;
end

%% Step 2 run ANN
tic
netLocked1 = fireNeurons(FirstDisp,Forecast,marginCost,scaleCost,IC);
tsim(1,2) = toc;
netLocked = netLocked1;
netLocked(netLocked1>=0.7) = 1;
netLocked(netLocked1<0.7) = 0;
% if Si==1 && timersOn
%     disp(['time for training is ' num2str(toc)])
% % elseif Si==3
% %     disp(['time for network propogation is ' num2str(toc)])
% end
tic
netAllLocked = ones(size(Locked));
netAllLocked(1,gens) = IC(gens)>0;
netAllLocked(2:end,gens) = netLocked';
QP = updateMatrices(Plant.OpMatB, IC, Time, scaleCost, marginCost,Forecast,renewgen,[]);
[GenDisp,~,Feasible] = DispatchQP(QP,netAllLocked);
if Feasible ~= 1
    netLocked(netLocked1>=0.8) = 1;
    netLocked(netLocked1<0.8) = 0;
    netAllLocked(2:end,gens) = netLocked';
    QP = updateMatrices(Plant.OpMatB, IC, Time, scaleCost, marginCost,Forecast,renewgen,[]);
    [GenDisp,~,Feasible] = DispatchQP(QP,netAllLocked);
%     [netAllDisp,~,Feasible] = DispatchQP(QPall,Organize,netAllLocked);
end
if Feasible ~= 1
    disp('no feasible solution from ANN')
end
% record.ANNtime((Si-2)*nS+1) = record.ANNtime((Si-2)*nS+1) + toc;
tsim(1,3) = toc;