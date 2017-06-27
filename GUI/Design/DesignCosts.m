function [MonthlyCosts,NetPresentCost] = DesignCosts(k)
global testSystems
GenDisp = testSystems(k).Design.GeneratorState;
Timestamp = testSystems(k).Design.Timestamp;
GenDisp = GenDisp(Timestamp>0,:);
Timestamp = Timestamp(Timestamp>0);
FuelInput = zeros(size(GenDisp));
MonthlyCosts = zeros(12,5); %12 months, each specifying: Electric Utility Demand charges, Electric Utility Use Charges, Fuel Costs, O&M, and Financing costs
nG = length(testSystems(k).Generator);
electricUtility = [];
for i = 1:1:nG
    skip = false;
    if strcmp(testSystems(k).Generator(i).Type,'Electric Generator') || strcmp(testSystems(k).Generator(i).Type,'CHP Generator')
        eff = testSystems(k).Generator(i).Output.Electricity;
        cap = testSystems(k).Generator(i).Output.Capacity*testSystems(k).Generator(i).Size;
    elseif strcmp(testSystems(k).Generator(i).Type,'Heater')
        eff = testSystems(k).Generator(i).Output.Heat; 
        cap = testSystems(k).Generator(i).Output.Capacity*testSystems(k).Generator(i).Size;
    elseif strcmp(testSystems(k).Generator(i).Type,'Utility') && strcmp(testSystems(k).Generator(i).Source,'Electricity')
        ElectricUse = GenDisp(:,i);
        electricUtility = i;
        skip = true;
    else skip = true;
    end
    if ~skip %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
        genEff = interp1(cap,eff,GenDisp(:,i));
        FuelInput(:,i) = GenDisp(:,i)./genEff;
        FuelInput(GenDisp(:,i)==0,i) = 0;
    end
end

scaleCost = updateGeneratorCost(Timestamp)*testSystems(k).optimoptions.Resolution;%% All costs were assumed to be 1 when building matrices
if ~isempty(electricUtility) && testSystems(k).Generator(electricUtility).VariableStruct.MinImportThresh<0
    if testSystems(k).Generator(electricUtility).VariableStruct.SellBackRate>0
        SellBack = testSystems(k).Generator(electricUtility).VariableStruct.SellBackRate*testSystems(k).optimoptions.Resolution;
    else
        SellBack = testSystems(k).Generator(electricUtility).VariableStruct.SellBackPerc*scaleCost(:,electricUtility);
    end
else SellBack = 0; %no sellback allowed
end

% startupcost = zeros(nnz(Timestamp),nG);
% for i = 1:1:nG
%     if isfield(testSystems(k).Generator(i).VariableStruct, 'StartCost')
%         nStartups = (GenDisp(2:end,i)>0).*(GenDisp(1:end-1,i)==0);
%         startupcost(2:end,i) = nStartups*testSystems(k).Generator(i).VariableStruct.StartCost;
%     end
% end

FuelCost = sum(FuelInput.*scaleCost,2);
ElecUseCost = max(0,ElectricUse).*scaleCost(:,electricUtility);
ElecSellBack = -min(0,ElectricUse).*SellBack;

n = length(testSystems(k).Costs.Equipment);
OandM = zeros(n,1);
FinanceMonthly = zeros(n,1);
for i = 1:1:n
    OandM(i) = testSystems(k).Costs.Equipment(i).OandM; %yearly maintenance cost
    LoanAmount = testSystems(k).Costs.Equipment(i).Cost*testSystems(k).Costs.Equipment(i).Financed/100; %principle loan amount
    MonthlyRate = testSystems(k).Costs.Equipment(i).LoanRate/12/100; 
    LoanTerm = testSystems(k).Costs.Equipment(i).LoanTerm;
    FinanceMonthly(i)=(MonthlyRate*LoanAmount/(1-(1+MonthlyRate)^(-LoanTerm*12))); % c*L/(1-(1+c)^(-n)) is the same as L*(c*(1+c)^n)/((1+c)^n -1)
    inflationPriceFactor(i,:)=(1.02.^(0:1:LoanTerm)); %annual inflaction factor for things like O &M, fuel and electric
end


%break into months and scale by % of month that was simulated
[y, m, ~] = datevec(Timestamp);
mUnique = unique(m); %months that were tested
for i = 1:1:length(mUnique)
    index = nonzeros((1:1:length(m))'.*(m==mUnique(i)));
    days = datenum([y(index(1)) mUnique(i)+1 1]) - datenum([y(index(1)) mUnique(i) 1]);
    fracOfMonth = length(index)/(days*24/testSystems(k).optimoptions.Resolution);
    MonthlyCosts(mUnique(i),1) = sum(FinanceMonthly);
    MonthlyCosts(mUnique(i),2) = sum(OandM)/12;
    MonthlyCosts(mUnique(i),3) = DemandCharge(ElectricUse(index),Timestamp(index),testSystems(k).Generator(electricUtility).VariableStruct);
    MonthlyCosts(mUnique(i),4) = sum(ElecUseCost(index))/fracOfMonth - sum(ElecSellBack(index))/fracOfMonth;
    MonthlyCosts(mUnique(i),5) = sum(FuelCost(index))/fracOfMonth;
end
%find net present cost
AnnualCost = sum(sum(MonthlyCosts))*12/length(mUnique);
NetPresentCost = NPC(1.02,AnnualCost*ones(20,1));%20 year net present cost without gas or electric escalators

function Cost = DemandCharge(Use,Timestamp,utility)
%find the highest demand charge of the month
day = weekday(Timestamp);
[Y,~,~,~,~,~] = datevec(Timestamp(1));
SumStart = datenum([Y(1),utility.SumStartMonth,utility.SumStartDay]);
WinStart = datenum([Y(1),utility.WinStartMonth,utility.WinStartDay]);
if Timestamp(1)>SumStart
    SumStart = datenum([Y(1)+1,utility.SumStartMonth,utility.SumStartDay]); %add a year
end
if Timestamp(1)>WinStart
    WinStart = datenum([Y(1)+1,utility.WinStartMonth,utility.WinStartDay]);% add a year
end
if WinStart<SumStart
    season = 'Sum';
else season = 'Win';
end
Rates = utility.(strcat(season,'Rates'))(:,2);
[~,~,~,H,~,~] = datevec(Timestamp);
OnPeak = zeros(length(day),1);
MidPeak = zeros(length(day),1);
OffPeak = zeros(length(day),1);
for t = 1:1:length(H) 
    j = utility.(strcat(season,'RateTable'))(day(t),H(t)+1);
    if j==1
        OnPeak(t) = Use(t);
    elseif j==2
        MidPeak(t) = Use(t);
    else
        OffPeak(t) = Use(t);
    end
end
Cost = max(OnPeak)*Rates(1);
Cost = Cost + max(MidPeak)*Rates(2);
Cost = Cost + max(OffPeak)*Rates(3);

function value=NPC(irr,cashflows)
value=0;
for i=1:length(cashflows)
    value=value+cashflows(i)/(irr)^i;
end