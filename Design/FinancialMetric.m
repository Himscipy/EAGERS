function result = FinancialMetric(heat, cool, elec, metric, params)
%FINANCIALMETRIC Calculate financial metric for given heating, cooling,
%electric demands.
% result = FINANCIALMETRIC(heat, cool, elec, metric, params)

bene = params.benefits;
cost = params.costs;
rate = params.rate;
time = 0:1:length(cost)-1;

switch lower(metric)
    case 'irr'
        tolerance = params.tolerance;
        npv = tolerance + 1;
        lowerBound = 0;
        upperBound = 1;
        growthStep = upperBound/2;
        rGuess = growthStep;
        boundedMode = false;
        while abs(npv) > tolerance
            npv = calc_npv(bene, cost, rGuess, time);
            if boundedMode
                % bounded search mode
                if npv > 0
                    lowerBound = rGuess;
                    rGuess = rGuess + (upperBound-rGuess)/2;
                else
                    upperBound = rGuess;
                    rGuess = rGuess + (lowerBound-rGuess)/2;
                end
            else
                % unbounded search mode
                if npv > 0
                    lowerBound = rGuess;
                    rGuess = upperBound;
                    upperBound = upperBound + growthStep;
                else
                    boundedMode = 1;
                    upperBound = rGuess;
                    rGuess = rGuess + (lowerBound-rGuess)/2;
                end
            end
        end
        result = rGuess;
    case 'lcoe'
        result = sum((Ct+Mt)/(1+rate)^time) / sum(Qt/(1+rate)^time);
    case 'npv'
        result = sum((Bt-Ct)/(1+rate)^time);
    case 'pbp'
        tolerance = params.tolerance;
        npv = tolerance + 1;
        lowerBound = 0;
        upperBound = 2^8;
        growthStep = upperBound/2;
        tGuess = 0:1:growthStep-1;
        boundedMode = false;
        while abs(npv) > tolerance
            npv = calc_npv(bene, cost, rate, tGuess);
            if boundedMode
                % bounded search mode
                if npv < 0
                    lowerBound = tGuess;
                    tGuess = tGuess + (upperBound-tGuess)/2;
                else
                    upperBound = tGuess;
                    tGuess = tGuess + (lowerBound-tGuess)/2;
                end
            else
                % unbounded search mode
                if npv < 0
                    lowerBound = tGuess;
                    tGuess = upperBound;
                    upperBound = upperBound + growthStep;
                else
                    boundedMode = 1;
                    upperBound = tGuess;
                    tGuess = tGuess + (lowerBound-tGuess)/2;
                end
            end
        end
        result = tGuess;
end

end


%% External functions
% Calculate net present value
function npv = calc_npv(b, c, r, t)
npv = sum((b-c) ./ (1+r) .^ t);
end
