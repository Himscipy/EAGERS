function varargout = AmbientTemperature(varargin)
%Provides current Ambient Temperature given a time of day (hours 1-24)
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'TemperatureData')
    SimSettings.TemperatureTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    SimSettings.TemperatureData = [22.7 22.3 22.0 21.7 21.4 22.0 22.7 23.3 24.6 25.9 27.3 27.8 28.4 29.0 29.1 29.2 29.3 28.4 27.4 26.5 25.5 24.6 23.7 23.2 22.7];
end
varargout{1} = interp1(SimSettings.TemperatureTimestamp,SimSettings.TemperatureData,mod(t/3600,24));
if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'TemperatureTimestamp';'TemperatureData'};
    varargout{2} = {'Time (in hours)','Temperature (C)'};
    varargout{3} = {'[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]','[22.7 22.3 22.0 21.7 21.4 22.0 22.7 23.3 24.6 25.9 27.3 27.8 28.4 29.0 29.1 29.2 29.3 28.4 27.4 26.5 25.5 24.6 23.7 23.2 22.7]'};
end