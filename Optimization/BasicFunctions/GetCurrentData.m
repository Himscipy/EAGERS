function Data = GetCurrentData(Date)
global  TestData
F = fieldnames(TestData.RealTimeData);
F = F(~strcmp('Timestamp',F));
nS = length(Date);
index = zeros(nS,1);
for t = 1:1:nS
    index(t) = nnz(TestData.RealTimeData.Timestamp<Date(t)+1e-6); %+1e-6 to avoid rounding problems
end
Data.Timestamp = Date;
for j = 1:1:length(F)
    if isstruct(TestData.RealTimeData.(F{j}))
        S = fieldnames(TestData.RealTimeData.(F{j}));
        for i = 1:1:length(S)
            Data.(F{j}).(S{i}) = TestData.RealTimeData.(F{j}).(S{i})(index,:);
        end
    else
        Data.(F{j}) = TestData.RealTimeData.(F{j})(index,:);
    end
end      