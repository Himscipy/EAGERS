function varargout = AmbientHumidity(varargin)
%Provides current Ambient Humidity given a time of day (seconds)
%Average hourly ambient humidity in July in DC
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'HumidityData')
    SimSettings.HumidityTimestamp = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    SimSettings.HumidityData = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];
end
varargout{1} = interp1(SimSettings.HumidityTimestamp,SimSettings.HumidityData,mod(t/3600,24));
if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'HumidityTimestamp';'HumidityData'};
    varargout{2} = {'Time (in hours)','Humidity (mass fraction)'};
    varargout{3} = {'[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]','[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]'};
end
end%Ends function AmbientHumidity