function Data = GetHistoricalData(Date)
global TestData
x1 = max(1,nnz(TestData.Timestamp<Date(1))); 
Data.Timestamp = Date;
F = fieldnames(TestData);
F = F(~strcmp('Timestamp',F));
F = F(~strcmp('HistProf',F));
if length(Date) == 1
    r = (Date - TestData.Timestamp(x1))/(TestData.Timestamp(x1+1) - TestData.Timestamp(x1));
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            if strcmp(F{j},'Hydro')
                S = S(~strcmp('Nodes',S));
            end
            for i = 1:1:length(S)
                Data.(F{j}).(S{i}) = (1-r)*TestData.(F{j}).(S{i})(x1,:) + r*TestData.(F{j}).(S{i})(x1+1,:);
            end
        elseif ~isempty(TestData.(F{j}))
            Data.(F{j}) = (1-r)*TestData.(F{j})(x1) + r*TestData.(F{j})(x1+1);
        end
    end
else
    x2 = max(1,nnz(TestData.Timestamp<Date(end))); 
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            if strcmp(F{j},'Hydro')
                S = S(~strcmp('Nodes',S));
            end
            for i = 1:1:length(S)
                Data.(F{j}).(S{i}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j}).(S{i})(x1:x2+1,:),Date);
            end
        elseif ~isempty(TestData.(F{j}))
            Data.(F{j}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j})(x1:x2+1),Date);
        end
    end
end