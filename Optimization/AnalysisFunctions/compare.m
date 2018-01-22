function compare
%% comparison for paper
MIPlant = load('Campus2_MixedInteger_Final.mat');
cQPlant = load('Campus_Non_MixedInteger_FullYear.mat');
%cQPlant = load('CampusRan.mat');
global Plant Si
Si = 2;
Plant = MIPlant.Plant;
endtime = 8713;
MIDisp = MIPlant.Plant.Dispatch.GeneratorState(2:endtime,:);
cQPDisp = cQPlant.Plant.Dispatch.GeneratorState(2:endtime,:);
ICMI = MIPlant.Plant.Dispatch.GeneratorState(1,:);
ICcQP = cQPlant.Plant.Dispatch.GeneratorState(1,:);


%find difference in dispatch setpoints
dispAble = [1,3,4,5,6,7,9,10,11,12,13,14,15,16,17];
ndisp = length(dispAble);
samesame = (MIDisp<cQPDisp+1).*(MIDisp>cQPDisp-1);
percSame = nnz(samesame(:,dispAble))/((endtime-1)*ndisp);

%find differences in dispatch on/off
gensOnOff = [3,4,5,6,9,10,11,12,13,17];%units that have non-zero lb
ngens = length(gensOnOff);
MIonoff = (MIDisp(:,gensOnOff)>0);
cQPonoff = (cQPDisp(:,gensOnOff)>0);

sameonoff = (MIonoff==cQPonoff);
percOnOff = nnz(sameonoff)/((endtime-1)*ngens);

%plot dispatches
% for i = 1:1:length(dispAble)
%     figure(i)
%     plot(MIDisp(:,dispAble(i)));
%     hold on
%     plot(cQPDisp(:,dispAble(i)));
%     title(strcat('disp for gen ',PlantcQP.Plant.Generator(dispAble(i)).Name));
%     legend('MI','cQP')
%     hold off
% end

%find difference in cost
DemandE = MIPlant.Plant.Dispatch.Demand.E(2:endtime);
Timestamp = MIPlant.Plant.Dispatch.Timestamp(2:endtime);
method = 'Dispatch';
[MICost,MImissedDem] = NetCostCalc(MIDisp,Timestamp,method);
[cQPCost,cQPmissedDem] = NetCostCalc(cQPDisp,Timestamp,method);
%cost difference
costDiff = cQPCost-MICost;
muCost = mean(costDiff);
sigmaCost = std(costDiff);
%individual method costs
mucQPCost = mean(cQPCost);
sigmacQPCost = std(cQPCost);
muMICost = mean(MICost);
sigmaMICost = std(MICost);
%missed demand due to difference between linear and non-linear chiller
%efficiency
muMIDem = mean(MImissedDem);
sigmaMIDem = std(MImissedDem);
mucQPDem = mean(cQPmissedDem);
sigmacQPDem = std(cQPmissedDem);

%plot 
%create gaussian of costs
nsigmas = 4;
dsigma = 1/8;
nbins = nsigmas*2*1/dsigma+1;

cQPcostBins = zeros(nbins*2-2,1);%bin by sigma/2, go out to six sigma +1 in each direction
MIcostBins = zeros(nbins*2-2,1);
%find bottom of bins
lbBincQP = mucQPCost-sigmacQPCost*dsigma/2-sigmacQPCost*nsigmas+[0:1:nbins-1]*sigmacQPCost*dsigma;
lbBinMI = muMICost-sigmaMICost*dsigma/2-sigmaMICost*nsigmas+[0:1:nbins-1]*sigmaMICost*dsigma;
%find top of bins
ubBincQP = mucQPCost+sigmacQPCost*dsigma/2-sigmacQPCost*nsigmas+[0:1:nbins-1]*sigmacQPCost*dsigma;
ubBinMI = muMICost+sigmaMICost*dsigma/2-sigmaMICost*nsigmas+[0:1:nbins-1]*sigmaMICost*dsigma;
for i = 1:1:nbins%count number in bins
    %find number in bins
    cQPcostBins((i-1)*2+1:i*2) = nnz(and((cQPCost>lbBincQP(i)),(cQPCost<=ubBincQP(i))));
    MIcostBins((i-1)*2+1:i*2) = nnz(and((MICost>lbBinMI(i)),(MICost<=ubBinMI(i))));
    %lbBincQP = lbBincQP+sigmacQPCost/2;%move up one sigma
    %ubBincQP = ubBincQP+sigmacQPCost/2;
end
%plot gaussian of cost
figure(1)
costscQP = sort([lbBincQP,ubBincQP]);
costsMI = sort([lbBinMI,ubBinMI]);
plot(costsMI, MIcostBins);
hold on
plot(costscQP,cQPcostBins);
title('Cost of Distribution of Dispatches')
xlabel('Cost in $ at each hour')
ylabel('number of hours in a year')
legend('mcQP', 'cQP');
xlim([0,2500]);
hold off    

%plot cost of onen day
iday = 4345;
figure(2)
plot(MICost(iday:iday+24))
hold on
plot(cQPCost(iday:iday+24))
hold off
title('Hourly Cost of Actual Dispatch for July 1st')
xlabel('Time in hours')
ylabel('Cost of disptach in $')
legend('mcQP', 'cQP')
JulyMICost = sum(MICost(iday:iday+24))
JulycQPCost = sum(cQPCost(iday:iday+24))
%show that different costs are related to different on/off, but that total
%cost is similar
Time = [1:24];
%plotComparison(MIDisp(iday:iday+24,:),Plant.Dispatch.Demand,Time,3)
%plotComparison(cQPDisp(iday:iday+24,:),Plant.Dispatch.Demand,Time,6)


%find out when fuel cells are at less than 99.5% of max power
FC1dwnMI = nnz(MIDisp(:,5)<1990);
FC2dwnMI = nnz(MIDisp(:,6)<1990);
FC1dwncQP = nnz(cQPDisp(:,5)<1990);
FC2dwncQP = nnz(cQPDisp(:,6)<1990);
FConoff = nnz(MIonoff(:,[3,4])==cQPonoff(:,[3,4]));
FCsame = nnz((MIDisp(:,[5,6])<cQPDisp(:,[5,6])+1).*(MIDisp(:,[5,6])>cQPDisp(:,[5,6])-1));
%percents
FConPerc = FConoff/(endtime*2);
FCsamePerc = FCsame/(endtime*2);

%find distance between generators of like types
%gas turbines 




%% compare concurrent run
%since the non-concurrent version shows different peaks, need to analyze
%horizon dispatches
%plot predicted dispatch costs
horizonCMI = MIPlant.Plant.Predicted.Cost(1:endtime-1);
horizonCcQP = MIPlant.Plant.Predicted.CostcQP(1:endtime-1);
horizonCcQP_noIC = zeros(endtime-1,1);
for i = 1:1:endtime-1
    horizonCcQP_noIC(i) = sum(NetCostCalc(cQPlant.Plant.Predicted.GenDisp(:,:,i),cQPlant.Plant.Predicted.Timestamp(:,i),method));
end
figure(9)
plot(horizonCMI)
hold on
plot(horizonCcQP)
hold off
title('Cost of Horizon Dispatch')
ylabel('Cost of Dispatch for 24 hour Horizon in $')
xlabel('Hour of the year')
xlim([0,8714])
legend('mcQP', 'cQP')
muhorizonCMI = mean(horizonCMI);
sigmahCMI = std(horizonCMI);
muhorizonCcQP = mean(horizonCcQP);
sigmahCcQP = std(horizonCcQP);

muCcQP_CMI = mean(horizonCMI-horizonCcQP);
muCcQPnoIC_CMI = mean(horizonCMI-horizonCcQP_noIC);

%compare starts
MIstarts = MIonoff(2:end,:)>MIonoff(1:end-1,:);
cQPstarts = cQPonoff(2:end,:)>cQPonoff(1:end-1,:);
nMIstarts = nnz(MIstarts);
ncQPstarts = nnz(cQPstarts);
end % End compare

function plotComparison(GeneratorDispatch,Demand,Time,fig)
%% plot dispatch
global Plant Si
%can't use Si here, because each step migh not be a whole number. Instead,
%round Si to the nearest whole number

dt = Time' - [0, Time(1:end-1)]';
Time = Time+(Si-2)*dt(1);
horizon = Plant.optimoptions.Horizon;
if horizon>24
    axStep = floor(horizon/10);
elseif horizon>12
    axStep=2;
else axStep=1;
end
axTick = ceil(Time(1)-dt(1)):axStep:Time(end);
axIndex = axTick;
while axIndex(end)>24
    axIndex(axIndex>24) = axIndex(axIndex>24)-24;
end

if isfield(Demand,'T')
    Demand = rmfield(Demand,'T');
end
S = {'E','H','C'};%Plant.optimoptions.Outputs;
for q = 1:1:length(S)
    name = {};
%     name2 = {};
    stor =[];
    plotBars =[];
    negBars =[];
    StoragePower =[];
    StorageState = [];
    FCcount = [];
    nLegend = zeros(1,length(Plant.Generator));
    for i = 1:1:length(Plant.Generator)
        if (strcmp(Plant.Generator(i).Type,'Thermal Storage') || strcmp(Plant.Generator(i).Type,'Electric Storage')) && isfield(Plant.Generator(i).QPform.output,S{q}) % energy storage
            StoragePower(:,end+1) = (GeneratorDispatch(1:end-1,i)-GeneratorDispatch(2:end,i))./dt;
            plotBars(:,end+1) = [0;max(0,StoragePower(:,end));];
            negBars(:,end+1) = -[0;max(0,-StoragePower(:,end));];
            name(end+1) = cellstr(Plant.Generator(i).Name);
            nLegend(i) = length(plotBars(1,:));
            stor(end+1) = i;
            %stor(end+1) = length(plotBars(1,:));
%             name(end+1) = {['Discharge ',Plant.Generator(i).Name]};
%             name2(end+1) = {['Charge ',Plant.Generator(i).Name]};
            if isfield(Plant.Generator(i).QPform.output,'E')
                bat = Plant.Generator(i).VariableStruct;
                DischEff = (bat.Voltage-bat.DischResist/10)/bat.Voltage;
                MaxDODcapacity = Plant.Generator(i).Size*(1-bat.MaxDOD/100);
                StorageState(:,end+1) = GeneratorDispatch(:,i)/DischEff+MaxDODcapacity;
            else
                StorageState(:,end+1) = GeneratorDispatch(:,i)/Plant.Generator(i).VariableStruct.DischargeEff;
            end
        elseif (strcmp(Plant.Generator(i).Type,'Utility')||~isempty(strfind(Plant.Generator(i).Type,'District'))) && isfield(Plant.Generator(i).QPform.output,S{q}) %utilities
            plotBars(:,end+1) = max(0,GeneratorDispatch(:,i));
            name(end+1) = cellstr(Plant.Generator(i).Name);
            if min(GeneratorDispatch(:,i))<0
                negBars(:,end+1)= -min(0,GeneratorDispatch(:,i));
                nLegend(i) = length(plotBars(1,:));
%                 name2(end+1) = cellstr('Sellback');
            end
        elseif isfield(Plant.Generator(i).QPform.output,S{q}) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S{q},'E')) %generators & renewables
            plotBars(:,end+1) = GeneratorDispatch(:,i)*Plant.Generator(i).QPform.output.(S{q})(1); %accounts for Hratio or such when gen has multiple outputs
            compName = cellstr(Plant.Generator(i).Name);
            name(end+1) = compName;
            if strcmp(compName{1}(1:2),'Fu') %if its a fuel cell
                FCcount(end+1) = i;
                nLegend(i) = length(plotBars(1,:));
            end
        end
    end
    %put storage at the end, this section should only exist for plotting
    %figures for comparison in a paper
    if isempty(stor)
        plotBars(:,end+1) = zeros(length(plotBars(:,1)),1);
    elseif stor(end)<length(plotBars(1,:))
        stori = nLegend(stor(1));
        plotBars = [plotBars(:,1:stori-1),plotBars(:,(stori+1):end),plotBars(:,stori)];
        name = {name{1:stori-1},name{stori+1:end},name{stori}};
        nLegend(stor) = nLegend(stor)+length(plotBars(1,(stori+1):end));
        nLegend = [nLegend(1:stor(1)-1),nLegend(stor(end)+1:end),nLegend(stor)];
    end
    %reorder fuel cells 1 and 2 to come first because they are almost
    %always on
    if ~isempty(FCcount)
        FC1 = nLegend(FCcount(1));
        FCend = nLegend(FCcount(end));
        name = {name{FC1:FCend},name{1:FC1-1},name{FCend+1:end}};
        plotBars = [plotBars(:,FC1:FCend),plotBars(:,1:FC1-1),plotBars(:,FCend+1:end)];
        nLegend(FCcount) = 0;
    end
    nLegend = nonzeros(nLegend);
    %% Plot
    if ~isempty(plotBars)
%         clf(fig+q-1)
        figure(fig+q-1)
        hold off
        plotTime = sort([Time(1)-1,Time,Time-.0001,Time(end)+dt(end)-.0001]);%make sure that there is a step wise characteristic of the plot
        posBars = zeros(length(plotBars(:,1))*2,length(plotBars(1,:)));
        posBars(1:2:end) = plotBars;
        posBars(2:2:end) = plotBars;
        h1 = area(plotTime,posBars,'Linestyle','none');
        if strcmp(S{q},'C')
            colormap cool
        elseif strcmp(S{q},'H')
            colormap autumn
        else
            colormap parula
        end
        legend(name,'Fontsize',16)
        hold on
        if ~isempty(negBars)
            OoM = log10(max(sum(plotBars,2)+sum(negBars,2)));
        else OoM = log10(max(sum(plotBars,2)));
        end
        if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
            Yspace = 10^(OoM-1);
            Ymax = 10^OoM;
        elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
            Yspace = 10^floor(OoM);
            Ymax = 10^ceil(OoM);
        elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
            Yspace = .5*10^floor(OoM);
            Ymax = .5*10^ceil(OoM);
        else  %count in increments of 2, 20, 200 or 2000 etc
            Yspace = .2*10^floor(OoM);
            Ymax = .2*10^ceil(OoM);
        end
        negTicks = 0;
        if ~isempty(negBars)
            negBarsPlot = zeros(length(negBars(:,1))*2,length(negBars(1,:)));
            negBarsPlot(1:2:end) = negBars;
            negBarsPlot(2:2:end) = negBars;
            h = area(plotTime,negBarsPlot,'Linestyle','none');%(2:end,:)'stacked','barwidth',1);
            negTicks = floor(min(min(negBars))/Yspace);
            if abs(negTicks)>3
                Yspace = 2*Yspace;
                negTicks = floor(min(min(negBars))/Yspace);
            end
            storageColors = {[.6 0 .6],[1 .8 .2],[.2 1 1]};%purple [.4 0 .6] if using jet, [.6 0 .6] if default map, orange, cyan
            for i = 1:1:length(nLegend)
                set(h1(nLegend(i)),'FaceColor',storageColors{i-(rem(i,3)-i)})
                set(h(i),'FaceColor',storageColors{i-(rem(i,3)-i)})
            end
        end
        Ymin = Yspace*negTicks;
        if isempty(Ymin)
            Ymin = 0;
        end
        pTicks = Ymax/Yspace;
        legend(name,'Fontsize',16,'Orientation','Horizontal','Location','NorthOutside','Box','off')
%         legend([name name2],'FontSize',16)
        xlim([Time(1)-dt(1), Time(end)])%[0,25]
        ylim([Ymin,Ymax])
        set(gca,'YTick',Ymin:Yspace:Ymax,'FontSize',14)
        set(gca, 'XTick',axTick,'XTickLabel', {axIndex})%{xticks(1:nS+1)})
        %% Storage
        if ~isempty(stor)
            [AX, H1, H2] =plotyy(0,0,[0,Time],StorageState');
            set(AX,{'ycolor'},{'k';'k'})
            set(H2(end),'Color','k','LineStyle','-','LineWidth',2)
            set(H2(1),'Color','k','LineStyle','-','LineWidth',2)
            xlim(AX(1),[Time(1)-dt(1),Time(end)])
            xlim(AX(2),[Time(1)-dt(1),Time(end)])
            set(AX(1),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',16)
            set(AX(2),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',16)
            ylim(AX(1),[Ymin,Ymax])
            set(AX(1),'YTick', Ymin:Yspace:Ymax)
            
            OoM = log10(max(max(StorageState)));
            if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
                Ymax = 10^OoM;
            elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
                Ymax = 10^ceil(OoM);
            elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
                Ymax = .5*10^ceil(OoM);
            else  %count in increments of 2, 20, 200 or 2000 etc
                Ymax = .2*10^ceil(OoM);
            end
            Yspace = Ymax/pTicks;
            Ymin = Yspace*negTicks;            
            ylim(AX(2),[Ymin,Ymax])
            set(AX(2),'YTick', Ymin:Yspace:Ymax)
            ylabel(AX(2),'State of Charge (kWh)','Color','k','FontSize',18)
        end
        xlabel('Time (hour)','Color','k','FontSize',18)
        ylabel('Generation (kW)','Color','k','FontSize',18)
        plotBars=[];
    end
end
end %End plotComparison