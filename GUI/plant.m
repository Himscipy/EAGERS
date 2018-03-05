classdef plant
    %both planning and dispatch tools use plant objects named Plant
    properties
        Name = 'New Plant';
        optimoptions = []; %options for optimization times/CHP/demands etc.
        Generator = [];
        Network = [];
    end
    methods 
        function obj = plant(name)
            obj.Name = name;
%             obj.Data.PlanningResult = [];
            %add optimization options
            obj.optimoptions.Interval = 31;
            obj.optimoptions.Horizon = 24;
            obj.optimoptions.Resolution = 1;
            obj.optimoptions.tspacing = 'constant';
            obj.optimoptions.excessHeat = true;
            obj.optimoptions.excessCool = false;
            obj.optimoptions.method = 'Dispatch';
            obj.optimoptions.MixedInteger = true;
            obj.optimoptions.SpinReserve = false;
            obj.optimoptions.SpinReservePerc = 0;
            obj.optimoptions.forecast = 'Perfect';
            obj.optimoptions.solver = 'quadprog';
            obj.optimoptions.mode = 'virtual';
            
%             obj.optimoptions.nsSmooth = 0;
%             obj.optimoptions.scaletime = 1;
%             obj.optimoptions.Topt = 3600;
%             obj.optimoptions.Tmpc = 600;
%             obj.optimoptions.thresholdSteps = 6;
            
            %add generators
            obj.Generator = [];
            obj.Generator.Type = [];
            obj.Generator.Name = [];
            obj.Generator.Source = [];
            obj.Generator.Output = [];
            obj.Generator.Size = [];
            obj.Generator.Enabled = [];
            obj.Generator.VariableStruct = [];
            
            %add Network structure
            obj.Network.name = '';
            obj.Network.Equipment = {};
            obj.Network.Location = [];
            %rest are networks
        end
    end
end
        