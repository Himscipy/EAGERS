function varargout = InternalHeatGains(varargin)
%Internal Gains per m^2 for lighting, solar insolation, and plug loads. Default is based off e+ Large Office
%% ** Note** the output needs to be multiplied by the area in m^2 to get the total W of internal heat gains
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'LightSched')
    SimSettings.InternalGainTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    SimSettings.LightSched = [.04619 .05 .05 .05 .05 .05 .1 .09238 .27714 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .4619 .27714 .27714 .18476 .18476 .09238 .04619];
    SimSettings.LightLoad = 9.6875;
    SimSettings.PlugSched = [.4 .4 .4 .4 .4 .4 .4 .4 .4 .9 .9 .9 .9 .8 .9 .9 .9 .9 .5 .4 .4 .4 .4 .4 .4];
    SimSettings.PlugLoad = 8.07007;
end
h = mod(t/3600,24);
varargout{1} = interp1(SimSettings.InternalGainTimestamp,SimSettings.LightSched*SimSettings.LightLoad,h) + SolarLoad + interp1(SimSettings.InternalGainTimestamp,SimSettings.PlugSched*SimSettings.PlugLoad,h);
if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'InternalGainTimestamp';'LightSched';'LightLoad';'PlugSched';'PlugLoad';};
    varargout{2} = {'Time (in hours)';'Lighting Schedule (fraction of peak)';'Lighting Peak load (W/m^2)';'Plug Load Schedule (fraction of peak plug load)';'Peak plug load (W/m^2)';};
    varargout{3} = {'[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]';...
                    '[.04619 .05 .05 .05 .05 .05 .1 .09238 .27714 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .83142 .4619 .27714 .27714 .18476 .18476 .09238 .04619]';...
                    '9.6875';...
                    '[.4 .4 .4 .4 .4 .4 .4 .4 .4 .9 .9 .9 .9 .8 .9 .9 .9 .9 .5 .4 .4 .4 .4 .4 .4]'...
                    '8.07007'};
end