function varargout = Occupancy(varargin)
%Occupancy per m^2  Default is based off e+ Large Office
%% ** Note** the output needs to be multiplied by the area in m^2 to get the total W of internal heat gains
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'Occupancy')
    SimSettings.OccupancyTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    SimSettings.Occupancy = [0.05 0 0 0 0 0 0 .1 .2 .95 .95 .95 .95 .5 .95 .95 .95 .95 .3 .1 .1 .1 .1 .05 .05];
    SimSettings.MaxOccupancy = 5.38195521e-2;
end
h = mod(t/3600,24);
varargout{1} = interp1(SimSettings.OccupancyTimestamp,SimSettings.Occupancy*SimSettings.MaxOccupancy,h);
if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'OccupancyTimestamp';'Occupancy';'MaxOccupancy';};
    varargout{2} = {'Time (in hours)';'Lighting Schedule (fraction of peak)';'Lighting Peak load (W/m^2)';'Plug Load Schedule (fraction of peak plug load)';'Peak plug load (W/m^2)';};
    varargout{3} = {'[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]';...
                    '[0.05 0 0 0 0 0 0 .1 .2 .95 .95 .95 .95 .5 .95 .95 .95 .95 .3 .1 .1 .1 .1 .05 .05]';...
                    '5.38195521e-2';};
end
end%Ends function Occupancy