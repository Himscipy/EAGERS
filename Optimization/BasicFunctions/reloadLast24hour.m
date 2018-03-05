function reloadLast24hour(Date,Resolution)
global TestData Last24hour
Last24hour = [];%re-load the previous 24 hours
Date = linspace(Date-1,Date,ceil(24/Resolution)+1)';
Last24hour.Timestamp = Date;
if Date(1)<TestData.Timestamp(1)
    Date = Date + ceil(TestData.Timestamp(1)-Date(1));%add a whole # of days
end
x1 = nnz(TestData.Timestamp<=Date(1)); 
F = fieldnames(TestData);
F = F(~strcmp('Timestamp',F));
F = F(~strcmp('RealTimeData',F));
x2 = max(1,nnz(TestData.Timestamp<Date(end))); 
for j = 1:1:length(F)
    if isstruct(TestData.(F{j}))
        S = fieldnames(TestData.(F{j}));
        for i = 1:1:length(S)
            if isnumeric(TestData.(F{j}).(S{i}))
                Last24hour.(F{j}).(S{i}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j}).(S{i})(x1:x2+1,:),Date);
            end
        end
    elseif ~isempty(TestData.(F{j})) && isnumeric(TestData.(F{j}))
        Last24hour.(F{j}) = interp1(TestData.Timestamp(x1:x2+1),TestData.(F{j})(x1:x2+1),Date);
    end
end