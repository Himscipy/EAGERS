function varargout = PowerDemandLookup(varargin)
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'PowerTime')
    SimSettings.PowerTime = [0 24];
    SimSettings.PowerDemand = [100 100]/100*SimSettings.NominalPower;
end
varargout{1} = interp1(SimSettings.PowerTime,SimSettings.PowerDemand,t/3600);

if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'PowerTime';'PowerDemand'};
    varargout{2} = {'Time (in hours)','Demand (% of Nominal)'};
    varargout{3} = {'[0 2 4 18 20 24]','[100 100 50 50 100 100]/100*SimSettings.NominalPower'};
end
end%Ends function PowerDemandLookup