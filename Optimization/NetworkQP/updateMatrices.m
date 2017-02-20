 function [QPall,Demand] = updateMatrices(QPall,Organize,IC,Time,scaleCost,marginCost,EC)
% % IC is the intial condition
% % fit refers to quadratic cost fit A or B
% % Locked is a matrix nSx nG with zeros when a generator is locked off
% % Stor is a variable to include storage in the optimization or not.

global Model_dir Plant DateSim UB Time
%function [QPall,Demand] = updateMatrices(QPall,Organize,IC,Time,scaleCost,marginCost,EC)
% IC is the intial condition
% fit refers to quadratic cost fit A or B
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Stor is a variable to include storage in the optimization or not.

ic = Plant.nodalVar.ic;
nodes = Plant.nodalVar.nodes; 
genNames = Plant.nodalVar.genNames;
nG = length(Plant.Generator);
nS = length(Time);

 Demand = updateForecast(DateSim,Time);%% function that creates demand vector with time intervals coresponding to those selected
 Demand = AccountForSelfDischarge(Demand,Time);

Outs = fieldnames(Plant.Data.Demand);    
for seq = 1:1:length(Outs)
    QP = QPall.(Outs{seq});% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb)
    thisSeq = Organize.(Outs{seq}).thisSeq;
    stor = Organize.(Outs{seq}).stor;
    utility = Organize.(Outs{seq}).utility;
    allStor = [stor];
    allGen = [thisSeq];
    allUtility = [utility];

    %updates demands: both beq & b (heating);
     for p = 1 %Placeholder to prevent running Heating
%    for p = 1:1:length(Plant.optimoptions.Outputs) %reserve nS steps for each output type
        mat = Organize.(Plant.optimoptions.Outputs{p}).Demand{1};
        index = Organize.(Plant.optimoptions.Outputs{p}).Demand{2};
        if strcmp(mat,'b')
            %QP.b(index) = -Demand.(Plant.optimoptions.Outputs{p});%Ax<beq, so b and A must be negative so you are producing more than enough
        else
            for t = 1:1:nS 
                for i = 1:1:nodes
                    if ~isempty(Plant.Network(i).demand)
                        QP.beq(((t-1)*nodes) + i + ic,1) = Plant.Data.Demand(Plant.Network(i).demand).E(t);
                    end 
                end
            end
        end  
    end


   %update initial condition
     for i = 1:1:nodes
        gen = Plant.Network(i).gen;
        for j = 1:1:length(gen)
            strf = strfind(gen{j},'.');
            genName = gen{j}(strf+1:end);
            I = find(strcmp(genName,genNames),1,'first');
             if ismember(I,[allGen, allStor]) && Organize.IC(I)~=0
                 QP.beq(Organize.IC(I)) = IC(I); 
                 QP.ub(Organize.States{I}(1)) = IC(I);
            end
            if ismember(I, allUtility)
                k = Organize.States{I}';
                if nnz(isinf(QP.ub(k)))>0
                    p = fieldnames(Plant.Generator(I).OpMatA.output); %utilities with infinite upper bound
                    QP.ub(k) = 100*max(sum(UB(~isinf(UB))),max(Demand.(p{1})));
                end
            end
        end
     end

    %update costs
    H = diag(QP.H);
    for i = 1:1:nodes
        gen = Plant.Network(i).gen;
        for j = 1:1:length(gen)
            s = strfind(gen{j},'.');
            genName = gen{j}(s+1:end);
            I = find(strcmp(genName,genNames),1,'first');
            if ismember(I,[allGen, allUtility])
                scaleC = [];
                k = Organize.States{I}';
                if ~isempty(k)
                    if isfield(Plant.Generator(I).OpMatA,'Ramp') %has ramping constraint, thus it has an initial condition in order to enforce this
                        k = k(2:end);%%first state corresponds to IC
                    end
                    for j = 1:1:floor(length(k)/nS)           
                        scaleC(end+1:end+nS,1) = scaleCost(:,I);   
                    end
                     H(k) = H(k).*scaleC;
                     QP.f(k) = QP.f(k).*scaleC;
                end
            elseif ismember(I,allStor)
                %% update storage costs
                %Storage States %SOC(t+1), charging power, upper buffer, lower buffer, self discharge 
                nX = Organize.States{I}; %first state corresponds to IC
                nSi = (length(nX)-length(Plant.Generator(I).OpMatA.states)+1); %final state of charge
                StorSize = Plant.Generator(I).OpMatA.X.ub;
                BuffSize = Plant.Generator(I).OpMatA.W.ub;
                if length(nX) == (4*nS+1) %If storage has 4 states plus ic
                    XlowBuff = 4:4:4*nS; %lower bound (W state) state #
                    XhighBuff = 5:4:1+(4*nS); %lower bound (W state) state #
%Done beleive we are using these currently
%                     IneqLowBuff = nS+1:2*nS; %lower bound (W state) inequality row
%                     IneqHighBuff = 2*nS+1:3*nS; %lower bound (W state) inequality row
                elseif length(nX) == (3*nS+1)
                    XlowBuff = 3:3:3*nS; %assuming no chargin penalty
                    XhighBuff = 4:3:1+(4*nS);
%                     IneqLowBuff = 1:nS;
%                     IneqHighBuff = nS+1:2*nS;
                     
                end   
                if ~isempty(EC)
                    %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
                    Out = fieldnames(marginCost);
                    for g = 1:1:length(Out)
                        if isfield(Plant.Generator(I).OpMatA,'Stor') && isfield(Plant.Generator(I).OpMatA.output,Out{g})
                            %% may need to increase the scale of the quadratic penalty to keep closer to original soln
                            %penalize more than 10% deviation in power output
                            PeakChargePower = Plant.Generator(I).OpMatA.Ramp.b(1);
                            dSOC_10perc = .1*PeakChargePower*Time(end); %energy (kWh) if charging at 10%
                            H(nX(nSi)) = -2*marginCost.(Out{g}).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                            QP.f(nX(nSi)) = -marginCost.(Out{g}).Min;%linear final value term loaded into SOC(t=nS)
    %                         QP.lb(nX(1:nSi)) = -EC(I);%change lb so that SOC = 0 coresponds to EC
    %                         QP.ub(nX(1:nSi)) = StorSize - EC(I);%change ub so that SOC = 0 coresponds to EC
                            %need to find the inequality rows and adjust b term to account for this offset
                            r = Organize.Inequalities{I};
    %                         QP.b(r(IneqLowBuff)) = -BuffSize + EC(I);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
    %                         QP.b(r(IneqHighBuff)) = StorSize-BuffSize - EC(I);%change lb so that SOC = 0 coresponds to EC
                            %need to change IC to IC-EC
    %                         QP.beq(Organize.IC(I)) = IC(I)-EC(I);
                        end
                    end
                else % dispatch optimization with final state & buffer states
                    Out = fieldnames(Plant.Generator(I).OpMatA.output);
                    if strcmp(Out{1},'H')
                        Max = 0.8*marginCost.H.Max;
                        Min = 0;
                    elseif strcmp(Out{1},'E')
                        Max = 1.25*marginCost.E.Max;
                        Min = 0.95*marginCost.E.Min;
                    elseif strcmp(Out{1},'C')
                        Max = 1.25*marginCost.C.Max;
                        Min = 0.65*marginCost.C.Min;
                    end
                    a1 = -Max; % fitting C = a1*SOC + a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                    a2 = (Max - Min)/(2*StorSize);
                    H(nX(nSi)) = 2*a2;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                    QP.f(nX(nSi)) = a1;%linear final value term loaded into SOC(t=nS)
                    if BuffSize>0
                        QP.f(nX(XlowBuff)) = Min;%this is the linear buffer term loaded into the upper buffer
                        QP.f(nX(XhighBuff)) = Min;%this is the linear buffer term loaded into the lower buffer
                        H(nX(XlowBuff)) = 2*(2*Max-Min)/(2*BuffSize);%this is the quadratic buffer term loaded into the upper buffer
                        H(nX(XhighBuff)) = 2*(2*Max-Min)/(2*BuffSize);%this is the quadratic buffer term loaded into the lower buffer    
                    end
                end
            end
        end
    end
    QP.H = diag(H);
    QPall.(Outs{seq})= QP;% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb) 
end 