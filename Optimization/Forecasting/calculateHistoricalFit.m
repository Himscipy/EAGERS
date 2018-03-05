global Plant
if isfield(Plant.Data,'Demand') && ~isempty(Plant.Data.Demand)
    F = fieldnames(Plant.Data.Demand);
    for i = 1:1:length(F)
        dem = F{i};
        h=waitbar(0,'Recalculating surface fit');

        Steps = round(1/(Plant.Data.Timestamp(2)-Plant.Data.Timestamp(1)));%points in 1 day of data
        Resolution =24/Steps;
        [~,m,~,hour,minutes,~] = datevec(Plant.Data.Timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
        m = [m(1); m; m(end)];%make m the same length as data
        hour = [max(0,floor(hour(1)-Resolution)); hour; min(24,hour(end)+floor((minutes(end)+60*Resolution)/60))];%make hour the same length as data
        minutes = [max(0,minutes(1)-60*Resolution); minutes; mod(minutes(end)+60*Resolution,60)];%make the same length as data
        months = sort(unique(m));

        monthNames = cellstr({'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';});
        BusDay = true(length(Plant.Data.Timestamp),1);
        daysAfter1_1_17 = floor(Plant.Data.Timestamp) - datenum([2017,1,1]);
        wd = 1 + mod(daysAfter1_1_17,7); %day of the week, sunday is 1, saturday is 7
        BusDay(wd==1) = false;
        BusDay(wd==7) = false;

        if length(months)<=4 %make 1 surface fit
            for k = 1:1:length(Plant.Data.Demand.(dem)(1,:))
                Z = Plant.Data.Demand.(dem)(:,k);
                Y = Plant.Data.Weather.Tdb;
                X = hour+minutes/60;
                if Steps>100 % too many data points for a surface fit, average each hour
                    r = ceil(Steps/24);
                    n2 = floor(length(X)/r);
                    X2 = zeros(n2,1);
                    Y2 = zeros(n2,1);
                    Z2 = zeros(n2,1);
                    for i=1:1:n2%turn while loop into for loop to make loop faster
                        X2(i)=sum(X((i-1)*r+1:i*r))/r;
                        Y2(i)=sum(Y((i-1)*r+1:i*r))/r;
                        Z2(i)=sum(Z((i-1)*r+1:i*r))/r;
                    end
                    X = X2;
                    Y = Y2;
                    Z = Z2;
                end
                Plant.Data.HistProf.(dem)(k).Annual = fit([X Y],Z,'lowess');
            end
        else %split surface fit by month and weekday/weekend
            Y = Plant.Data.Weather.Tdb;
            X = hour+minutes/60;
            for k = 1:1:length(Plant.Data.Demand.(dem)(1,:))
                Z = Plant.Data.Demand.(dem)(:,k);
                for i = 1:1:length(months)
                    waitbar((i-1)/length(months),h,strcat('Calculating surface fit for ',dem,' for the month of ',monthNames{months(i)}));
                    if Steps>150 % too many data points for a surface fit, average each 10 min
                        r = ceil(Steps/(24*6));
                        n2 = floor(length(X)/r);
                        X2=zeros(n2,1);
                        Y2=zeros(n2,1);
                        Z2=zeros(n2,1);
                        for j=1:1:n2%turn the while loop into a for loop to make loop faster
                            X2(j)=sum(X((j-1)*r+1:j*r))/r;
                            Y2(j)=sum(Y((j-1)*r+1:j*r))/r;
                            Z2(j)=sum(Z((j-1)*r+1:j*r))/r;
                        end
                        X = X2;
                        Y = Y2;
                        Z = Z2;
                    end
                    % business days
                    H1 = (m==months(i)) & BusDay; %points in the correct month
                    Z1 = Z(H1);%remove the days that are non business days
                    Y1 = Y(H1);
                    X1 = X(H1);%cannot use nonzeros function because there are some elements in X1 that are 0
                    % weekends/holidays
                    H2 = (m==months(i)) & ~BusDay;
                    Z2 = Z(H2);
                    Y2 = Y(H2);
                    X2 = X(H2);
                    Plant.Data.HistProf.(dem)(k).(strcat(monthNames{months(i)},'WeekDay')) = fit([X1, Y1],Z1,'lowess');%'poly23');%'loess');
                    Plant.Data.HistProf.(dem)(k).(strcat(monthNames{months(i)},'WeekEnd')) = fit([X2, Y2],Z2,'lowess');%'poly23');%'loess');
                end
            end
        end
        close(h)
    end
end