function varargout = PowerDemandLookup(varargin)
global SimSettings modelParam
t = varargin{1};
if ~isfield(SimSettings,'PowerTime')
    SimSettings.PowerTime = [0 24];
    SimSettings.PowerDemand = [100 100];
end
varargout{1} = interp1(SimSettings.PowerTime,SimSettings.PowerDemand,t/3600)/100*modelParam.NominalPower;

if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'PowerTime';'PowerDemand'};
    varargout{2} = {'Time (in hours)','Demand (% of Nominal)'};
    varargout{3} = {'[0 4 8 24]','[100 100 50 50]'};
end