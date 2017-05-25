function [fcast,err] = SNIWPEForecast(histData, histTime, hzn, ...
    prevFcast, prevErr, beta)
%[fcast,err] = SNIWPEForecast(histData, histTime, hzn, prevFcast,
%       prevErr, beta)
%Seasonal Naive Incorporating Weighted Progressive Error forecast
%   INPUTS:
%       histData            - list of historical data values (most recent
%                               value at last index)
%       histTime            - list of hours histData is sampled at
%       hzn                 - list of horizons to forecast
%       prevFcast           - forecast at the last timestep (may be input
%                               as NaN if no previous forecast exists)
%       prevErr             - previous error value used
%       beta                - value for determining bias correction weights
%   OUTPUTS:
%       fcast               - list of forecasted values
%       err                 - error used in correction
%   
%   SNIWPE forecasts are based on past data using the equation:
%       f_t+h|t = y_(t-s+h) - e_t * w
%   f_t+h|t is the forecast h timesteps into the future given data up
%       to time t.
%   y_(t-s+h) is the data point one period (or season) behind the data
%       point to be forecasted (a period could be between one hour and one
%       year).
%   e_t is the value of the dynamic error at time t:
%       e_t = f_t|(t-1) - y_t + e_(t-1)
%   w is a vector of weights determined by beta. For linear weighting,
%       w = max(1-slp.*((1:1:p)-1), 0)

%% Built-in function parameters
p = 24; % seasonality (or period) of data (in hours)

%% Check inputs
assert(beta >= 0 && beta <= 1, 'beta should be between 0 and 1')
[s1,s2] = size(histData);
if s2 ~= 1
    assert(s1==1, 'histData should be a vector; not a matrix')
    histData = histData';
end
[s1,s2] = size(histTime);
if s2 ~= 1
    assert(s1==1, 'histTime should be a vector; not a matrix')
    histTime = histTime';
end
[s1,s2] = size(hzn);
if s2 ~= 1
    assert(s1==1, 'hzn should be a vector; not a matrix')
    hzn = hzn';
end

%% Determine weighting
hitZeroVal = max(beta*p, histTime(2));
slp = -1 / (hitZeroVal-histTime(1));
w = max(1+slp.*(hzn-hzn(1)), 0);

%% Forecast
if ~isnan(prevFcast)
    err = prevFcast - histData(end) + prevErr;
else
    err = 0;
end
nDays = hzn(end)/24; % number days in hzn
eaDay = length(hzn) / nDays; % number indices in each day
naiveFcast = zeros(length(hzn),1);
% Use: interp1(x, values, xquery)
for i = 1:1:ceil(nDays)
    rng = 1+floor((i-1)*eaDay):1:min(floor(i*eaDay), length(hzn)); % range
    xNaive = mod(hzn(rng), 24);
    yNaive = interp1(histTime,histData,xNaive);
    naiveFcast(rng) = yNaive;
end

fcast = naiveFcast - err * w;

end