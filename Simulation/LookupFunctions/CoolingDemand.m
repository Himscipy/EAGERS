function varargout = CoolingDemand(varargin)
global SimSettings
t = varargin{1};
if ~isfield(SimSettings,'CoolingDemand')
    SimSettings.CoolingTime = [0 24];
    SimSettings.CoolingDemand = [100 100]/100*SimSettings.NominalCooling;
end
h = mod(t/3600,24);
varargout{1} = interp1(SimSettings.CoolingTime,SimSettings.CoolingDemand,h);
if length(varargin)==2
    %output what the schedule variables are
    varargout{1} = {'CoolingTime';'CoolingDemand'};
    varargout{2} = {'Time (in hours)','CoolingDemand (kW)'};
    varargout{3} = {'[0 4 8 24]','[100 100 50 50]/100*SimSettings.NominalCooling'};
end
end%Ends function CoolingDemand