%TEST_FINANCIALMETRIC

netCost = [
    95000   1       0
    3       15      2
    10      6       3
    200     1       4
    -5      1       4
    8       1       5
];
netBaselineCost = [
    700     1       3
];
rate = 0.01638;
lifetimeYrs = 20;
tolerance = 1e-2;
iterLimit = 100;
metric = 'irr';

params.cost = netCost;
params.baseCost = netBaselineCost;
params.rate = rate;
params.time = [lifetimeYrs 4];
params.tolerance = tolerance;
params.iterLimit = iterLimit;

%% Get result
result = FinancialMetric(metric, params);

%% Post-results actions and plot preparation
switch metric
    case 'irr'
        rate = result.irr;
        if result.flag
            resultStr = 'Did not converge';
        else
            resultStr = num2str(rate);
        end
        fprintf('IRR:\t%s\n', resultStr)
        % Iterations of Rate and NPV
        figure(1)
        plot(result.iterVecNpv, '-o')
        title('Iterations: Net Present Value')
        xlabel('iteration')
        ylabel('npv')
        figure(2)
        plot(result.iterVecR, '-o')
        title('Iterations: Discount Rate')
        xlabel('iteration')
        ylabel('rGuess')
        % NPV vs. Rate
        nQuery = 40;
        qStep = (result.rHighest - result.rLowest) / nQuery;
        queryR = result.rLowest : qStep : result.rHighest;
        queryNpv = zeros(size(queryR));
        for i = 1:1:length(queryR)
            queryNpv(i) = sum((result.bene-result.cost) ./ ...
                (1+queryR(i)) .^ result.time);
        end
    case 'npv'
        fprintf('NPV:\t%s\n', num2str(result.npv))
end

% Non-discounted and discounted cash flows
beneCash = zeros(length(result.time),1);
costCash = zeros(length(result.time),1);
beneNpv = zeros(length(result.time),1);
costNpv = zeros(length(result.time),1);
for i = 1:1:length(result.time)
    beneCash(i) = beneCash(max(i-1,1)) + ...
        result.bene(i);
    costCash(i) = costCash(max(i-1,1)) + ...
        result.cost(i);
    beneNpv(i) = beneNpv(max(i-1,1)) + ...
        result.bene(i)/(1+rate)^result.time(i);
    costNpv(i) = costNpv(max(i-1,1)) + ...
        result.cost(i)/(1+rate)^result.time(i);
end

%% Plot
% Non-discounted cash flow
figure(3)
plot(result.time, -beneCash, '-o')
hold on
plot(result.time, -costCash, '-o')
title('Net Cash Flow Along Project Lifetime')
legend({'baseline','scenario A'})
xlabel('Time (yrs)')
ylabel('Cash Flow')
hold off
% Discounted cash flow
figure(4)
plot(result.time, -beneNpv, '-o')
hold on
plot(result.time, -costNpv, '-o')
title('Net Present Value for Varying Project Lifetime')
legend({'baseline','scenario A'})
xlabel('Time (yrs)')
ylabel('NPV')
hold off
% NPV vs. Rate
figure(5)
plot(queryR, queryNpv, '-o')
title('Effect of Rate on NPV')
xlabel('Rate')
ylabel('NPV')
