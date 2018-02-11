function DesignCosts(k,Years,Costs)
global testSystems
GenDisp = testSystems(k).Design.GeneratorState;
FuelInput = zeros(size(GenDisp));

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
    else
        skip = true;%skip chillers whose input shows up as electric or heat generation, skip storage devices
    end
    if ~skip %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
        genEff = interp1(cap,eff,GenDisp(:,i));
        FuelInput(:,i) = GenDisp(:,i)./genEff;
        FuelInput(GenDisp(:,i)==0,i) = 0;
    end
end

scaleCost = updateGeneratorCost(testSystems(k).Design.Timestamp,testSystems(k).Generator)*testSystems(k).optimoptions.Resolution;%% All costs were assumed to be 1 when building matrices
FuelCost = sum(FuelInput.*scaleCost,2);
if ~isempty(electricUtility) && testSystems(k).Generator(electricUtility).VariableStruct.MinImportThresh<0
    if testSystems(k).Generator(electricUtility).VariableStruct.SellBackRate>0
        SellBack = testSystems(k).Generator(electricUtility).VariableStruct.SellBackRate*testSystems(k).optimoptions.Resolution;
    else
        SellBack = testSystems(k).Generator(electricUtility).VariableStruct.SellBackPerc*scaleCost(:,electricUtility);
    end
else
    SellBack = 0; %no sellback allowed
end
ElecUseCost = max(0,ElectricUse).*scaleCost(:,electricUtility);
ElecSellBack = -min(0,ElectricUse).*SellBack;

n = length(Costs);
OandM = zeros(n,1);
FinanceMonthly = zeros(n,1);
CapitalCost = zeros(n,1);
for i = 1:1:n
    if ~isempty(Costs(i).Name)
        CapitalCost(i) = Costs(i).Cost;
        OandM(i) = Costs(i).OandM; %yearly maintenance cost
        LoanAmount = Costs(i).Cost*Costs(i).Financed/100; %principle loan amount
        MonthlyRate = Costs(i).LoanRate/12/100; 
        LoanTerm = Costs(i).LoanTerm;
        FinanceMonthly(i)=(MonthlyRate*LoanAmount/(1-(1+MonthlyRate)^(-LoanTerm*12))); % c*L/(1-(1+c)^(-n)) is the same as L*(c*(1+c)^n)/((1+c)^n -1)
    end
end

MonthlyCosts = zeros(0,6); 
D = datevec(testSystems(k).Design.Timestamp(1));
m = 0;
fracOfMonth = (datenum([D(1) D(2)+1 1]) - testSystems(k).Design.Timestamp(1))/(datenum([D(1) D(2)+1 1])-datenum([D(1) D(2) 1]));
while datenum([D(1) D(2)+m 1])<testSystems(k).Design.Timestamp(end)
    index = max(1,nnz(testSystems(k).Design.Timestamp<=datenum([D(1) D(2)+m 1]))):1:nnz(testSystems(k).Design.Timestamp<datenum([D(1) D(2)+m+1 1]));
    m = m+1;
    startupcost = zeros(nG,1);
    for i = 1:1:nG
        if isfield(testSystems(k).Generator(i).VariableStruct, 'StartCost')
            nStartups = nnz((GenDisp(index(2:end),i)>0) & (GenDisp(index(1:end-1),i)==0));
            startupcost(i) = nStartups*testSystems(k).Generator(i).VariableStruct.StartCost;
        end
    end
    MonthlyCosts(end+1,1) = sum(FinanceMonthly);
    MonthlyCosts(m,2) = sum(OandM)/12;
    MonthlyCosts(m,3) = sum(startupcost);
    MonthlyCosts(m,4) = DemandCharge(ElectricUse(index),testSystems(k).Design.Timestamp(index),testSystems(k).Generator(electricUtility).VariableStruct);
    MonthlyCosts(m,5) = sum(ElecUseCost(index))/fracOfMonth - sum(ElecSellBack(index))/fracOfMonth;
    MonthlyCosts(m,6) = sum(FuelCost(index))/fracOfMonth;
    fracOfMonth = min(1,(testSystems(k).Design.Timestamp(end) - datenum([D(1) D(2)+m 1]))/(datenum([D(1) D(2)+m 1])-datenum([D(1) D(2)+m-1 1])));
end
testSystems(k).Costs.Design = MonthlyCosts;
%find net present cost
Months = 12*Years;
Escalator = ones(Months,6);
eEscalator = ProjectUtilityCost('electric','Commercial',Years,'WA')';
eEscalator = eEscalator/eEscalator(1);
gEscalator = ProjectUtilityCost('gas','Commercial',Years,'WA')';
gEscalator = gEscalator/gEscalator(1);
Escalator(:,4) = eEscalator;
Escalator(:,5) = eEscalator;
Escalator(:,6) = gEscalator;
mC = zeros(Months,1);
mC(1:m) = sum(MonthlyCosts(1:m,2:end),2);
mC(m+1:12) = sum(mC(1:m))/m;
for n = m+1:1:12%make sure there are at least 12 months
    MonthlyCosts(end+1,:) = sum(MonthlyCosts(1:m,:),1)/m;
end
for i = m+1:1:Months
    mC(i) = sum(MonthlyCosts(rem(i-1,12)+1,2:end).*Escalator(i,2:end));%monthly costs without financing charges
end
testSystems(k).Costs.NPC = NPC((1+testSystems(k).Costs.DiscountRate/100),mC+sum(FinanceMonthly));
mC(1) = mC(1)+sum(CapitalCost);%put all capital cost in 1st month
testSystems(k).Costs.ProjectedMonthlyCosts = mC;


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
else
    season = 'Win';
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
    value=value+cashflows(i)/(irr)^(i/12);
end
