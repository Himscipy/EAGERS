function Data = GetHistoricalData(Date)
global TestData
Data.Timestamp = Date;
if Date(1)<TestData.Timestamp(1)
    Date = Date + ceil(TestData.Timestamp(1)-Date(1));%add a whole # of days
end
x1 = nnz(TestData.Timestamp<=Date(1)); 
F = fieldnames(TestData);
F = F(~strcmp('Timestamp',F));
F = F(~strcmp('RealTimeData',F));
if length(Date) == 1
    r = (Date - TestData.Timestamp(x1))/(TestData.Timestamp(x1+1) - TestData.Timestamp(x1));
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                if isnumeric(TestData.(F{j}).(S{i}))
                    Data.(F{j}).(S{i}) = (1-r)*TestData.(F{j}).(S{i})(x1,:) + r*TestData.(F{j}).(S{i})(x1+1,:);
                end
            end
        elseif ~isempty(TestData.(F{j})) && isnumeric(TestData.(F{j}))
            Data.(F{j}) = (1-r)*TestData.(F{j})(x1) + r*TestData.(F{j})(x1+1);
        end
    end
else
    x2 = max(1,nnz(TestData.Timestamp<Date(end))); 
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                if isnumeric(TestData.(F{j}).(S{i}))
                    Data.(F{j}).(S{i}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j}).(S{i})(x1:x2+1,:),Date);
                end
            end
        elseif ~isempty(TestData.(F{j})) && isnumeric(TestData.(F{j}))
            Data.(F{j}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j})(x1:x2+1),Date);
        end
    end
end