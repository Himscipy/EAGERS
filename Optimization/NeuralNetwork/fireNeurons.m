  function LockedNet = fireNeurons(FirstDisp,GenDemand,marginCost,scaleCost,IC)
global Plant dischEff Si DateSim
global lockedNet lockederror net %these are created in this function and used in each subsequent step
global Last24hour
%input: inputs for neuron dispatch and network training
%output: dispatch from network
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = [Time - [0; Time(1:end-1)]]';
nG = length(Plant.Generator);
nS = length(Time);
allgens = 1:nG;
Out = Plant.optimoptions.Outputs;
renews = []; %number of renewables
gens = []; %generators
stor = []; %storage
UB = zeros(1,nG);
LB = zeros(1,nG);
QPsingle = Plant.OpMatB;
f = zeros(1,nG);%linear costs of each comp
H = zeros(1,nG);%quadratic costs of each comp
Hratio = zeros(1,nG);
Cratio = zeros(1,nG);
for i = 1:1:nG%find renewables and upper and lower bounds
    states = QPsingle.Organize.States{i};
    if strcmp(Plant.Generator(i).Source, 'Renewable')
        renews(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = 0;
    elseif ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
        stor(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = UB(i).*.2;
    elseif ~strcmp(Plant.Generator(i).Type,'Utility') && ~strcmp(Plant.Generator(i).Type,'Heater')%no unit commitment for utilities, storage, renewables, or heaters
        gens(end+1) = i;
        UB(i) = Plant.Generator(i).Size;
        LB(i) = min(0,sum(QPsingle.lb(states)));
        f(i) = Plant.OpMatB.f(states(1));%linear costs
        H(i) = Plant.OpMatB.H(states(2),states(2));%quadratic costs
        if isfield(Plant.Generator(i).Output, 'H')
            Hratio(i) = Plant.Generator(i).Output.H(end);
        end
    else
        UB(i) = inf;
        LB(i) = 0;
    end
end
gensandstor = [gens,stor];   
ngs = length(gensandstor);
renewgen = zeros(nS,nG);
ngens = length(gens);

if Si==1
    training = true;
else training = false;
end
%% create training data for 7 days at different start times
if training %train the first time through
    Tdays = min(14, floor(length(Plant.Data.Demand.E)/(24*4)));%number of training data sets (Training days). You can only train for the number of days you have data for
    DateSim1 = DateSim;
    histdataTime = datevec(Plant.Data.Timestamp(1));
    DateSimVec = datevec(DateSim);
    DateSim = datenum([histdataTime(1), DateSimVec(2:5), 0]);%make sure you are training over the same year you have data for, don't include seconds
    Last24hour1 = Last24hour;
    tic
    GenDisp = zeros(nS*Tdays,nG);%preallocate recording of solution
    scaleCostT = zeros(nS*Tdays,length(scaleCost(1,:)));%preallocate recording of cost
    FirstDispT = zeros(nS*Tdays,length(FirstDisp(1,:)));%preallocate recording of first dispatch (FitA)
    RenewE = zeros(nS*Tdays,1);%preallocate recording of renewable gen
    PredictDispatch = [];
    maxGenDem = 0;
    DemandT1day = GenDemand;
    for i = 1:1:length(Out)
        DemandT.(Out{i}) = zeros(nS*Tdays, length(GenDemand.(Out{i})(1,:)));%this is demand Total: all the demand including what is provided by renewables over all the training days
    end
    for day = 0:1:Tdays-1%run for number of training days
       %% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
       DateSim = DateSim+1;%Plant.Data.Timestamp(day+1);
        if isempty(PredictDispatch)
            PredictDispatch = ones(length(Time)+1,1)*IC;
        end
        scaleCost = updateGeneratorCost(Time/24 + DateSim); %% All feedstock costs were assumed to be 1 when building matrices 
        marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
        
        %the demand should be historical data + gaussian noise (-renewable
        %generation if it is added to the QP)
        for i =1:1:length(renews)%renewable generation
            renewgen(:,renews(i))= RenewableOutput(Plant.Generator(renews(i)).VariableStruct,DateSim,Time,'Predict');%Add forecast of renewable power @ resulution equal to the Time vector
        end
        %rgen = sum(renewgen,2);%this assumes that all renewable generation goes to electric demand
        for j = 1:1:length(Out)
            histData = Plant.Data.Demand.(Out{j})(int16((DateSim-Plant.Data.Timestamp(1))*24*4+[1:nS].*dt*4'));%assumes data is every 15 minutes
            sizeHistData = size(histData);
            if sizeHistData(1) == 1 %if format is in a single row, change it to a single collumn
                Last24hour.(Out{j}) = histData;
                histData = histData';
            else Last24hour.(Out{j}) = histData';
            end
            sigma = std(histData)/5;
            qpdemand = max(histData + randn(size(histData))*sqrt(sigma),0);%historical data plus gaussian noise
            DemandT.(Out{j})(1+nS*day:nS*(1+day),:) = qpdemand;
            DemandT1day.(Out{j}) = qpdemand;
        end
        QP = updateMatrices(Plant.OpMatB,IC,Time,scaleCost,marginCost,DemandT1day,renewgen,[]); %update fit B matrices
        
        %% Step 1 Determine initial dispatch
        Locked = true(nS+1,nG);
        for i = 1:1:nG
            if ~Plant.Generator(i).Enabled
                Locked(:,i) = 0;
            end
        end
        [FirstDisp,~,Feasible] = DispatchQP(QP,Locked);
        if ~(Feasible==1)%% hopefully not here
            disp('error in training: initial dispatch was infeasible. Defaulting to previous dispatch');
            FirstDisp= PredictDispatch;
        end
        %% Step 2 Dispatch step by step
%         if Plant.optimoptions.MixedInteger
%             tic
%             OptimalState = StepByStepDispatch(Forecast,Renewable,scaleCost,dt,IC,'initially constrained',FirstDisp);
%             tsim(1,2) = toc;
%         else
            OptimalState = FirstDisp;
%         end
        %% Start with optimal dispatch, and check if feasible
        tic
        marginCost = updateMarginalCost(OptimalState,scaleCost,Time);
        QP = updateMatrices(Plant.OpMatB,IC,Time,scaleCost,marginCost,DemandT1day,renewgen,[]); %update fit B matrices
        for i = 1:1:nG
            if QP.Organize.Dispatchable(i) ==1
                Locked(OptimalState(:,i)==0,i)=false;
            end
        end
%         if Plant.optimoptions.MixedInteger
%             [DayDisp, ~, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B
%             if Feasible ~=1
%                 [DayDisp, ~, Feasible] = FindFeasible(QP,Locked);
%                 if ~(Feasible==1)%% hopefully not here
%                     disp('error: Cannot Find Feasible Dispatch');
%                 end
%             end
%         else
            DayDisp = FilterGenerators(QP,OptimalState,Locked,[0;Time/24] + DateSim);
%         end
        %record
        GenDisp(1+nS*day:nS*(day+1),:) = DayDisp(2:end,:);
        scaleCostT(1+nS*day:nS*(day+1),:) = scaleCost;
        FirstDispT(1+nS*day:nS*(1+day),:) = FirstDisp(2:end,:);
        PredictDispatch = DayDisp;%[DayDisp(ceil(.3*nS):end,:);DayDisp(1:ceil(.3*nS)-1,:)];
        maxGenDem = max(max(GenDemand.E),maxGenDem);
        RenewE(1+nS*day:nS*(day+1),:) = sum(renewgen,2);
        Last24hour.Timestamp = DateSim+Time./24;
    end
    GenDisp = [GenDisp(1,:);GenDisp];
    scaleCost = scaleCostT;
    FirstDisp = [FirstDispT;FirstDispT(end,:)];
    for i = 1:1:length(Out)
        GenDemand.(Out{i}) = DemandT.(Out{i});
    end
    Time = repmat(Time,Tdays,1);
    dt = repmat(dt,1,Tdays);
    nS = length(Time);
    DateSim = DateSim1;
    Last24hour = Last24hour1;
end

%% create train and dispatch a locked network
lockInputs = zeros(nS, 8*ngens+1+length(stor)+length(Out));
for i = 1:1:length(Out)
    if strcmp(Out,'E')
        lockInputs(:,i) = (GenDemand.(Out{i})-RenewE)./(max(GenDemand.(Out{i})));
    else
        lockInputs(:,i) = GenDemand.(Out{i})./(max(GenDemand.(Out{i})));
    end
end
col = length(Out)+1;%collumn for input
lockInputs(:,col) = (scaleCost(:,1).*dt')./max(Plant.Generator(1).VariableStruct.SumRates(:,1)).*dt';%cost of utility electrical
lockInputs(:,col+1:col+length(stor)) = (FirstDisp(1:end-1,stor)-FirstDisp(2:end,stor))./(ones(nS,1)*max(FirstDisp(1:end-1,stor)-FirstDisp(2:end,stor)));%SOC of battery at previous step
col = col+length(stor);
for i = 1:1:ngens
    j = gens(i);
    j1 = QPsingle.Organize.States{i};
    if ~isempty(j1)
        if ~isfield(Plant.Generator(j).OpMatB.output,'C') || ~Plant.optimoptions.sequential %if its not a chiller or the chillers are run concurrently pull costs from OneStep.E
            lockInputs(:,col) = scaleCost(:,2).*f(j)./(max(scaleCost(:,2))*max(f)-min(scaleCost(:,2))*min(f));
            lockInputs(:,col+1) = scaleCost(:,2).*H(j)./(max(scaleCost(:,2))*max(H));%minH is 0
            lockInputs(:,col+2) = Hratio(j)./max(Hratio);
            col = col+3;
        else %if its a chiller have cost be proportional to the energy it uses up
            lockInputs(:,col) = scaleCost(:,1).*(Cratio(j))./(max(Cratio)*max(scaleCost(:,1)));
            col = col+1;
            %lockInputs(:,col+i) = (scaleCost(:,2).*Plant.OpMatB.f(j1(1)).*dt')./max(scaleCost(:,2).*Plant.OpMatB.f(j1(1)).*dt');%linear cost
            %lockInputs(:,col+ngens+i) = (scaleCost(:,2).*Plant.OpMatB.H(j1(2),j1(2)).*dt')./max(scaleCost(:,2).*Plant.OpMatB.H(j1(2),j1(2)).*dt');%quadratic cost
        end
%         if isfield(Plant.Generator(j).VariableStruct, 'StartCost') %heaters don't have this
%             lockInputs(:,col+2*ngens+i) = (scaleCost(:,2).*Plant.Generator(j).VariableStruct.StartCost.*dt')./max(scaleCost(:,2).*Plant.Generator(j).VariableStruct.StartCost.*dt');%onCost
%         end
        lockInputs(:,col) = min(FirstDisp(1:end-1,j)+Plant.Generator(j).OpMatB.Ramp.b(1),UB(j))./UB(j);%upper bound
        lockInputs(:,col+1) = max(FirstDisp(1:end-1,j)-Plant.Generator(j).OpMatB.Ramp.b(1),LB(j))./UB(j);%lower bound
        lockInputs(:,col+2) = FirstDisp(2:end,j)./UB(j);
        lockInputs(:,col+3) = (FirstDisp(1:end-1,j)-LB(j))./UB(j);
        lockInputs(:,col+4) = [FirstDisp(3:end,j);FirstDisp(1,j)]./UB(j);
        col = col+4;
    end
end

%prevent Nans
lockInputs(isnan(lockInputs)) = 0;
%removeCols = true(1,length(lockInputs(1,:))); %if no variance, remove that collumn
% for i = 1:1:length(lockInputs(1,:))
%     if nnz(lockInputs(:,i)~=0)==0
%         removeCols(i) = false;
%     end
% end
%lockInputs = lockInputs(:,removeCols);

if training %train it the first time through, use it later
    lockedNet = Neural_Network(length(lockInputs(1,:)),ngens,'classify',1);
    %train locked network
    [lockedNet, lockederror] = trainNetwork(lockedNet, (GenDisp(2:end,gens)>0)',lockInputs);
    nS = nS/Tdays;
end
%use the network to find which generators to lock off or on
LockedNet = forward(lockedNet, lockInputs(1:nS,:));
% LockedNet(LockedNet>=0.95) = 1;
% LockedNet(LockedNet<0.95) = 0;

%% create train and dispatch dispatch network
% geninputs = zeros(nS, length(lockInputs(1,:))+3*length(stor)+length(LockedNet(:,1)));
% geninputs(:,1:length(lockInputs(1,:))) = lockInputs;
% index = length(lockInputs(1,:));
% %index = index+ngens*2;
% %add the discharge efficiency, and marginal costs
% for i = 1:1:length(stor)
%     j = stor(i);
%     geninputs(:,index+i) = dischEff(j);
%     %geninputs(:,index+length(stor)+i) = FirstDisp(1:end-1,j);%SOC at previous step
%     geninputs(:,index+length(stor)*2+i) = (marginCost.E.Min)./max(max(scaleCost));%margin cost at time t
%     geninputs(:,index+length(stor)*3+i) = marginCost.E.Max./max(max(scaleCost));
% end
% index = index+length(stor)*3;
% 
% allLocked = LockedNet;
% %if a generator is locked off, change all inputs to 0
% geninputs(:,index+1:end) = LockedNet';
% % if sum(sum(abs(lockederror)))>nS/10 %if incorrect for more than 5%, only use as an input
% %     allLocked = ones(ngens,nS);
% % end
%     
% if training
%     if sum(sum(abs(lockederror)))>nS/10%if incorrect for more than 5%, use the training Locked config
%         allLocked = (GenDisp(2:end,gens)>0)';
%     end
%     %rescale the generator dispatch for the upper/lower bound
%     net = Neural_Network(length(geninputs(1,:)),ngens);%ngs);
%     for i = 1:1:ngens%ngs
%         j = gens(i);
%         geninputs1 = geninputs.*(ones(length(geninputs(1,:)),1)*allLocked(i,:))';
%         componentnet = Neural_Network(length(geninputs(1,:)),1);
%         GenDispTrain = max((GenDisp(2:end,j)'-LB(j))./(UB(j)-LB(j)),0);
%         [componentnet, sqrerror] = trainNetwork(componentnet, GenDispTrain, geninputs1);
%         if nnz(sqrerror>.005)>nS*ngs/8 %if more than an 8th the results have high error try training again
%             shuffle = round(rand(nS,1).*nS);
%             shuffle(shuffle==0) = 1;
%             geninputsshuff = geninputs1(shuffle, :);
%             GenDispTrainshuff = GenDispTrain(:,shuffle);
%             [componentnet2,sqrerror2] = trainNetwork(componentnet, GenDispTrainshuff, geninputsshuff);
%             if mean(mean(sqrerror2))<mean(mean(sqrerror))
%                 componentnet = componentnet2;
%                 sqrerror = sqrerror2;
%             end
%         end
%         net.Wlayer1(:,i) = componentnet.Wlayer1;
%     end
% end
%  
% if training && timersOn
%     disp(['time for training neural network is ' num2str(toc)])
% end
% 
% % tic
% netDisp = zeros(nS,ngs);
% for i = 1:1:ngens %add in storage later
%     j = gens(i);
%     geninputs1 = geninputs.*(ones(length(geninputs(1,:)),1)*allLocked(i,:))';
%     netDisp(:,i) = (geninputs1*net.Wlayer1(:,i)).*(UB(j)-LB(j))+LB(j).*allLocked(i,:)';
%     netDisp(netDisp(:,i)<LB(j)/2,i) = 0;
% end
% % if training && timersOn
% %     disp(['time for dispatching neural network is ' num2str(toc)])
% % end
% if training && nnz(abs(netDisp-GenDisp(2:end,gensandstor))>0.1)>0 %if there is a .1kW difference between the dispatch from the net and dispatch from StepByStep
%     disp('Neural Network dispatch differs from OprimalState. training needs revision.')
% end

