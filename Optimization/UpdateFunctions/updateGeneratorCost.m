function scaleCost = updateGeneratorCost(Timestamp,Gen)
nS = length(Timestamp);
Source = {zeros(length(Gen),1)};
for i = 1:1:length(Gen)
    if strcmp(Gen(i).Type,'Utility')
        Source(i) = {Gen(i).Source};
        utility = Gen(i).VariableStruct;
        if strcmp(Source(i), 'Electricity')
            day = weekday(Timestamp);
            [Y,M,D] = datevec(Timestamp(1));
            WinStart = datenum([Y,utility.WinStartMonth,utility.WinStartDay]);
            if M>utility.WinStartMonth || (M==utility.WinStartMonth && D>=utility.WinStartDay)
                SumStart = datenum([Y+1,utility.SumStartMonth,utility.SumStartDay]);
            else
                SumStart = datenum([Y,utility.SumStartMonth,utility.SumStartDay]);
            end
            Rate = zeros(nS,1);
            [~,~,~,H,~,~] = datevec(Timestamp);
            H(H==0) = 24; %zeroth hour is column 1 in rate matrix
            for t = 1:1:length(H) 
                if Timestamp(t)>=SumStart && Timestamp(t)<WinStart
                    Rate(t) = utility.SumRates(utility.SumRateTable(day(t),H(t),1));
                else
                    Rate(t) = utility.WinRates(utility.WinRateTable(day(t),H(t),1));
                end
            end
            Utility(i).Rate = Rate;%rate in $/kWh
        elseif strcmp(Source(i), 'NG')|| strcmp(Source(i), 'Diesel')
            %gas utility set up as daily prices for a leap year
            date1 = datevec(utility.Timestamp(1));
            utility.Timestamp(end+1) = utility.Timestamp(end)+1;
            utility.Rate(end+1) = utility.Rate(1);
            year1 = date1(1);
            datenow = datevec(Timestamp);
            interpDate = datenum([year1*ones(length(Timestamp),1), datenow(:,2:end)]);
            Utility(i).Rate = interp1(utility.Timestamp,utility.Rate,interpDate)/293.1; %interpolate & convert gas rate from $/MMBTU to $/kWh;
        else
            %% need to add something for district heating/cooling
        end
    end
end

scaleCost = zeros(nS,length(Gen));
for i = 1:1:length(Gen)
    if strcmp(Gen(i).Type,'Utility')
        scaleCost(:,i) = Utility(i).Rate;
    elseif isempty(strfind(Gen(i).Type,'Storage')) && ~strcmp(Gen(i).Source, 'Renewable') && ~strcmp(Gen(i).Type, 'Hydro') %not storage or renewable or hydro
        Uindex = find(strcmp(Source,Gen(i).Source),1);
        if isempty(Uindex)
            scaleCost(:,i) = 1; %no utility (don't scale costs)
        else
            scaleCost(:,i) = Utility(Uindex).Rate;
        end
    end
end