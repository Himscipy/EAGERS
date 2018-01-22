function errors = Compare2EplusCompHand(building, weather, date, varargin)
%COMPARE2EPLUSCOMPHAND Compare EAGERS and EnergyPlus results using
%ComparisonHandler object.
%   errors = Compare2EplusCompHand(building, weather, date[, flags, args])
%
%   FLAG OPTIONS:
%   'comp'      Comparative plot. Info for such plotting should be
%               immediately following 'comp' and should be a cell array
%               containing name of variable in building.VariableStruct to
%               be studied, list of values to study, and optionally a cell
%               array containing the names of plots to keep (default: keep
%               all plots). For this option, ERRORS will be output as a
%               matrix; rows corresponding to load type and columns
%               corresponding to different iterations.
%               Example:
%               compVar =       'Resistance';
%               compValues =    [700, 800, 900];
%               plotsToKeep =   {'Heating', 'Cooling'};
%               compInfo = {compVar, compValues, plotsToKeep};
%               errors = Compare2EplusCompHand(building, weather, date, ...
%               'comp', compInfo);
%
%   'stack'     Stacked bar plot. Each bar is represents the span of one
%               timestep and is split up by demand category.
%
%   See also COMPARISONHANDLER

%% Handle input DATE
[s1, ~] = size(date);
if s1 == 1
    date = date';
end

%% Use ComparisonHandler object
compHand = ComparisonHandler(building, weather, date);
errors = compHand.run(varargin{:});

end
