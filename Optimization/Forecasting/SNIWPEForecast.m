function Forecast = SNIWPEForecast(Date)
global Last24hour
Forecast.Timestamp = Date;
%% initialize
alpha = .85;
beta = 0.9; % Take this out after training is added!
Outs = fieldnames(Last24hour.Demand);
for S = 1:1:length(Outs)
    if ~isfield(Last24hour,'prevErr') || ~isfield(Last24hour.prevErr,Outs{S}) || isnan(Last24hour.prevFcast.(Outs{S}))
        Last24hour.prevErr.(Outs{S}) = 0;
        Last24hour.prevFcast.(Outs{S}) = Last24hour.Demand.(Outs{S});
    else
        Last24hour.prevErr.(Outs{S}) = Last24hour.prevFcast.(Outs{S}) - Last24hour.Demand.(Outs{S})(end) + Last24hour.prevErr.(Outs{S});
    end
end
Dvec = Last24hour.Timestamp;
nDays = ceil(Date(end)-Last24hour.Timestamp(end));
w = [linspace(alpha,0,nnz(Date<(Last24hour.Timestamp(end)+beta)))';zeros(nnz(Date>=(Last24hour.Timestamp(end)+beta)),1)];
n = length(Dvec)-1;
for d = 1:1:nDays
    Dvec(end+1:end+n) = Last24hour.Timestamp(2:end) + d;
end
for S = 1:1:length(Outs)
    data = Last24hour.Demand.(Outs{S});
    for d = 1:1:nDays
        data(end+1:end+n) = Last24hour.Demand.(Outs{S})(2:end);
    end
    Forecast.Demand.(Outs{S}) = interp1(Dvec,data,Date);
    Forecast.Demand.(Outs{S})(1:length(w)) = Forecast.Demand.(Outs{S})(1:length(w)) + w*Last24hour.prevErr.(Outs{S}); 
end