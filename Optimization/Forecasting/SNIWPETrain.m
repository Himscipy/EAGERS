function optBeta = SNIWPETrain(histData)

%% Get MAPE for widely spaced beta values
betas = (0 : 0.1 : 1)';


%% Hone in on the correct beta value with a smaller step size
betas = (bestBeta-0.1 : 0.01 : bestBeta+0.1)';

end
