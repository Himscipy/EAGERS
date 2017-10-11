%% comparison for paper
MIPlant = load('Campus_MixedInteger_7.mat');
cQPlant = load('CampusRan.mat');
global Plant Si
Si = 1;
Plant = MIPlant.Plant;
endtime = 700;
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
dsigma = 1/10;
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
legend('complementary QP', 'Modified cQP');
hold off    

%plot cost of first day
figure(2)
plot(MICost(1:24))
hold on
plot(cQPCost(1:24))
hold off
title('Hourly Cost of Actual Dispatch for January 1')
xlabel('time in hours')
ylabel('Cost of disptach in $')
legend('complementary QP', 'Modified cQP')
%show that different costs are related to different on/off, but that total
%cost is similar
Time = [1:24];
plotDispatch2(MIDisp(1:25,:),Plant.Dispatch.Demand,Time,3)
plotDispatch2(cQPDisp(1:25,:),Plant.Dispatch.Demand,Time,6)


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
figure(9)
plot(MIPlant.Plant.Predicted.CostcQP)
hold on
plot(MIPlant.Plant.Predicted.Cost)
hold off
title('Cost of Horizon Dispatch')
ylabel('Cost of Dispatch for 24 hour Horizon in $')
xlabel('Hour of the year')
legend('complementary QP', 'Modified cQP')



