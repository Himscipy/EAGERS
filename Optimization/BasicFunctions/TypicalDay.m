function Z = TypicalDay(h,Timestamp,Data)
Steps = round(1/(Timestamp(2)-Timestamp(1)));%points in 1 day of data
Resolution =24/Steps;
[~,m,~,hour,minutes] = datevec(Timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
m = [m(1); m; m(end)];%make m the same length as data
hour = [max(0,floor(hour(1)-Resolution)); hour; min(24,hour(end)+floor((minutes(end)+60*Resolution)/60))];%make hour the same length as data
months = sort(unique(m));
%% fit temperature data
Z = ones(12,24);
for i = 1:1:length(months)
    if ~isempty(h)
        waitbar(i/length(months),h,'Calculating Monthly Profiles');
    end
    Total = zeros(24,1);
    points = zeros(24,1);
    for j = 1:1:24
        X = (m==i) & (hour == j-1);
        Total(j) = sum(Data(X));
        points(j) = nnz(Data(X));
    end
    Z(months(i),:) = (Total./points)';
end
if length(months)<12
    if nnz(months==1)==1 && nnz(months==12)==1
        j =1;
        for i= 2:1:11
            if nnz(months==i)==1
                j=j+1;
            else Z(i,:) = Z(j,:);
            end
        end
    else
        for i= 1:1:min(months)-1
            Z(i,:) = Z(min(months),:);
        end
        for i= max(months)+1:1:12
            Z(i,:) = Z(max(months),:);
        end
    end
end
end%Ends function TypicalDay