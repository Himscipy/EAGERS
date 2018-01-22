function [Cost,Eimbalance,Himbalance] = NetCostCalc(Var1,Timestamp,method)
global Plant
letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';};
if strcmp(method,'Dispatch')
    Dispatch = Var1;
    Input = 0*Dispatch;
    Eimbalance = 0*Timestamp;
    Himbalance = 0*Timestamp;
    for i = 1:1:length(Plant.Generator)
        skip = false;
        if ~isempty(Plant.Generator(i).Output)
            cap = Plant.Generator(i).Output.Capacity*Plant.Generator(i).Size;
        end
        if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
            eff = Plant.Generator(i).Output.Electricity;
        elseif strcmp(Plant.Generator(i).Type,'Chiller') 
            eff = Plant.Generator(i).Output.Cooling;
            COP = interp1(cap,eff,Dispatch(:,i));
            skip = true;
            if strcmp(Plant.Generator(i).Source,'Electricity')
                E_use = Plant.Generator(i).QPform.constDemand.E*(Dispatch(:,i)>0);
                Output = Dispatch(:,i);
                for j = 1:1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end)))
                    E_use = E_use - min(Output,Plant.Generator(i).QPform.(letters{j}).ub(end))*Plant.Generator(i).QPform.output.E(j,end);
                    Output = max(0,Output - min(Output,Plant.Generator(i).QPform.(letters{j}).ub(end)));
                end
                Eimbalance = Eimbalance + (Dispatch(:,i)./COP - E_use);
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                H_use = Plant.Generator(i).QPform.constDemand.H*(Dispatch(:,i)>0);
                Output = Dispatch(:,i);
                for j = 1:1:nnz(~cellfun('isempty',Plant.Generator(i).QPform.states(:,end)))
                    H_use = H_use - min(Output,Plant.Generator(i).QPform.(letters{j}).ub(end))*Plant.Generator(i).QPform.output.H(j,end);
                    Output = max(0,Output - min(Output,Plant.Generator(i).QPform.(letters{j}).ub(end)));
                end
                Himbalance = Himbalance + (Dispatch(:,i)./COP - H_use);
            end
        elseif strcmp(Plant.Generator(i).Type,'Heater')
            eff = Plant.Generator(i).Output.Heat;  
        else skip = true;
        end
        if ~skip %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
            Input(:,i) = Dispatch(:,i)./interp1(cap,eff,Dispatch(:,i));
            if nnz(Dispatch(:,i))<length(Timestamp) %if you have steps where the dispatch is zero, prevent Nan
                OffSteps = (Dispatch(:,i)==0);
                Input(OffSteps,i) = 0;
            end
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            Input(:,i) = Dispatch(:,i);
        end
    end
    %% ___ %%% alternate: uses the cost in the optimization (FitB)
%     dt = Timestamp(2:end)-Timestamp(1:end-1);
%     Cost = zeros(length(Timestamp)-1,1);
%     StorPower =[];
%     for i = 1:1:length(Plant.Generator)
%         if nnz(Dispatch(:,i))>0
%             if ~isempty(Plant.Generator(i).QPform.states)
%                 states = Plant.Generator(i).QPform.states(:,end);
%             end
%             if isfield(Plant.Generator(i).QPform,'Stor')
%                 StorPower(:, end+1) = Dispatch(1:end-1,i) - Dispatch(2:end,i); %no cost
%             elseif isfield(Plant.Generator(i).QPform,'constCost') %all of these cost terms need to be scaled later on
%                 I = Plant.Generator(i).QPform.(states{1}).ub(end);
%                 b1 = Plant.Generator(i).QPform.(states{1}).f(end);
%                 b2 = Plant.Generator(i).QPform.(states{2}).f(end);
%                 b3 = Plant.Generator(i).QPform.(states{2}).H(end);
%                 b4 = Plant.Generator(i).QPform.constCost;
%                 Cost = Cost + (min(Dispatch(2:end,i),I)*b1 + max(Dispatch(2:end,i)-I,0)*b2 + max(Dispatch(2:end,i)-I,0).^2*b3+b4).*dt.*scaleCost(:,i);
%             elseif~isempty(states) %single state generators (linear cost term) & utilities
%                 b1 = Plant.Generator(i).QPform.(states{1}).f(end);
%                 Cost = Cost + Dispatch(2:end,i)*b1.*scaleCost(:,i).*dt;
%                 %need to handle case of sell-back
%             elseif ~isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
%                 %renewable
%             end
%         end
%     end  
elseif strcmp(method,'Input')
    Input = Var1;
    Dispatch = 0*Input;
    Dispatch(Input>1) = 1;
end
nS = nnz(Timestamp);
nG = length(Plant.Generator);
dt = (Timestamp(2:nS) - Timestamp(1:nS-1))*24;
scaleCost = updateGeneratorCost(Timestamp(1:nS-1)).*(dt*ones(1,nG));%% All costs were assumed to be 1 when building matrices
stor =[];
startupcost = zeros(nS-1,nG);
for i = 1:1:nG
    if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
        nStartups = (Dispatch(2:nS,i)>0).*(Dispatch(1:nS-1,i)==0);
        startupcost(:,i) = nStartups*Plant.Generator(i).VariableStruct.StartCost;
    end
    if isfield(Plant.Generator(i).QPform,'Stor')
       stor(end+1) = i; 
    end
end
scaleCost(:,stor)=0; %remove cost of storage
Cost = sum(Input(2:nS,1:nG).*scaleCost + startupcost,2);