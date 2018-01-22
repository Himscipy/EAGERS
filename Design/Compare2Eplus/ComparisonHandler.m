classdef ComparisonHandler < handle
% ComparisonHandler Handles comparison between EAGERS and EnergyPlus.
%   ComparisonHandler(building, weather, date, nFigs)
%   BUILDING is a struct containing building information.
%   WEATHER is a struct containing weather information.
%   DATE is a list of datenum values.
%
%   See also COMPARE2EPLUSCOMPHAND

%   To add a plot field:
%   1   constructor:    add to "obj.figTitles" list.
%   2   eagers_results: add to "loads" list.
%   3   eplus_results:  add to "loads" list.
    
    properties
        bldg
        wthr
        date
        mtr
        eagers
        ePlus
        figTitles
        nFigs
    end
    
    methods
        % Constructor
        function obj = ComparisonHandler(bldg, wthr, date)
            obj.bldg = bldg;
            obj.wthr = wthr;
            obj.date = date;
            obj.mtr = obj.load_meter();
            obj.eagers = obj.run_eagers();
            obj.ePlus = obj.load_eplus_results();
            obj.figTitles = {...
                'Equipment', ...
                'Interior Lighting', ...
                'Exterior Lighting', ...
                'Cooling', ...
                'Heating', ...
                'Fan Power', ...
                'HVAC', ...
                'Total'};
            obj.nFigs = length(obj.figTitles);
            obj.check_nfigs()
        end
        
        % Load EnergyPlus meter file
        function meter = load_meter(obj, mtrName)
            global Model_dir
            if nargin < 2
                mtrName = strcat(obj.bldg.Name, '.mat');
            end
            mtrFile = fullfile(Model_dir, 'Design', 'Compare2Eplus', mtrName);
            mtrLoaded = load(mtrFile);
            meter = mtrLoaded.mtr;
        end
        
        % Run EAGERS and store the results
        function loads = run_eagers(obj)
            [equipment, ...
                interiorLighting, ...
                exteriorLighting, ...
                cooling, ...
                heating, ...
                fanPower, ...
                otherLoads] = ...
                BuildingProfile(obj.bldg, obj.wthr, obj.date, [],obj.mtr);
            cooling = cooling / obj.bldg.VariableStruct.COP_C;
            heating = heating / obj.bldg.VariableStruct.COP_H;
            hvac = heating + cooling + fanPower;
            total = equipment + interiorLighting + exteriorLighting + ...
                hvac + otherLoads;
            loads = {...
                equipment, ...
                interiorLighting, ...
                exteriorLighting, ...
                cooling, ...
                heating, ...
                fanPower, ...
                hvac, ...
                total};
        end
        
        % Results of running EnergyPlus
        function loads = load_eplus_results(obj)
            loads = {...
                obj.mtr.InteriorEquipmentElectricity, ...
                obj.mtr.InteriorLightsElectricity, ...
                obj.mtr.ExteriorLightsElectricity, ...
                obj.mtr.CoolingElectricity, ...
                obj.mtr.HeatingElectricity, ...
                obj.mtr.FansElectricity, ...
                obj.mtr.ElectricityHVAC, ...
                obj.mtr.ElectricityNetFacility};
        end
        
        % Check number of figures reported
        function check_nfigs(obj)
            nfGood = obj.nFigs==length(obj.eagers) && ...
                obj.nFigs==length(obj.ePlus);
            errText = sprintf(strcat('The following should all have an', ...
                ' equal number of elements:\n', ...
                '  - figure titles\n', ...
                '  - EAGERS results\n', ...
                '  - EnergyPlus results'));
            assert(nfGood, errText)
        end
        
        % Main function (public)
        function errors = run(obj, varargin)
            errors = obj.elec_load_plot(varargin{:});
        end
        
        % Electric load plotting using flags
        function percErrors = elec_load_plot(obj, varargin)
            specialModes = {'comp', 'stack'};
            mode = 'default';
            for i = 1:1:length(specialModes)
                specialIndex = find(strcmp(varargin, specialModes{i}));
                if specialIndex
                    mode = specialModes{i};
                    info = varargin{specialIndex + 1};
                end
            end
            switch mode
                case 'comp'
                    percErrors = obj.plot_comparison(info);
                case 'stack'
                    percErrors = obj.plot_stack(info);
                otherwise
                    percErrors = obj.plot_default();
            end
            obj.order_figs()
        end
        
        % Plot in default mode
        function percErrors = plot_default(obj)
            [xE, dOfY, dt, dtE, xi, xf] = obj.time_constants();
            percErrors = zeros(obj.nFigs, 1);
            for iFig = 1:1:obj.nFigs
                figTitle = obj.figTitles{iFig};
                fig = figure(iFig);
                obj.add_plot_elems(fig, figTitle)
                obj.step_plot(xE, obj.ePlus{iFig}(xi:xf), 'b')
                hold on
                obj.step_plot(dOfY, obj.eagers{iFig}, 'r')
                hold off
                legend({'EnergyPlus', 'EAGERS'})
                PercError = 100 * (sum(obj.eagers{iFig}(2:end).*dt) - ...
                    sum(obj.ePlus{iFig}.*dtE)) / sum(obj.ePlus{iFig}.*dtE);
                fprintf('%s electric load:\t%f\n', figTitle, PercError)
                percErrors(iFig) = PercError;
            end
        end
        
        % Plot in comparison mode
        function percErrors = plot_comparison(obj, compInfo)
            compVar     = compInfo{1};  % character array
            compVal     = compInfo{2};  % number array
            plotsToKeep = {};           % character cell array
            if length(compInfo) > 2
                plotsToKeep = compInfo{3};
            end
            nCompVals = length(compVal);
            percErrors = zeros(obj.nFigs, nCompVals);
            for iVal = 1:1:nCompVals
                obj.display_iteration_text(compVar, compVal, iVal)
                obj.update_simulation_results(compVar, compVal, iVal)
                percErrors(:, iVal) = obj.comparison_iteration(obj.ePlus, ...
                    obj.eagers, plotsToKeep, iVal);
                fprintf('\n')
            end
        end
        
        % Display comparison plot iteration to the console
        function display_iteration_text(~, compVar, compVal, iVal)
            fprintf('-----Iteration %i-----\n', iVal)
            if iscell(compVal)
                fprintf('--> %s = %s\n', compVar, compVal{iVal})
            else
                fprintf('--> %s = %f\n', compVar, compVal(iVal))
            end
        end
        
        % Update simulation results for comparison mode
        function update_simulation_results(obj, compVar, compVal, iVal)
            if strcmp(compVar, 'mtr')
                obj.mtr = obj.load_meter(compVal{iVal});
                obj.ePlus = obj.load_eplus_results();
            else
                obj.bldg.VariableStruct.(compVar) = compVal(iVal);
                obj.eagers = obj.run_eagers();
            end
        end
        
        % One comparison plot iteration
        function percErrors = comparison_iteration(obj, ePlus, eagers, ...
                plotsToKeep, iVal)
            [xE, dOfY, dt, dtE, xi, xf] = obj.time_constants();
            percErrors = zeros(obj.nFigs, 1);
            shouldPlot = sum(strcmpi(plotsToKeep, 'none')) == 0;
            for iFig = 1:1:obj.nFigs
                figTitle = obj.figTitles{iFig};
                keepCondition = isempty(plotsToKeep) ...
                    || sum(strcmp(plotsToKeep, figTitle));
                if shouldPlot && keepCondition
                    fig = figure(iFig);
                    if iVal == 1
                        obj.add_plot_elems(fig, figTitle)
                        obj.step_plot(xE, ePlus{iFig}(xi:xf))
                        legend('EnergyPlus')
                    end
                    hold on
                    dispName = sprintf('EAGERS %i', iVal);
                    obj.step_plot(dOfY, eagers{iFig}, 'DisplayName', dispName)
                    hold off
                end
                pErr = 100 * (sum(eagers{iFig}(2:end).*dt) - ...
                    sum(ePlus{iFig}(xi:xf).*dtE)) / ...
                    sum(ePlus{iFig}(xi:xf).*dtE);
                fprintf('%s electric load:\t%f\n', figTitle, pErr)
                percErrors(iFig) = pErr;
            end
        end
        
        % Plot one comparison graph
        function plot_comparison_graph(obj)
            
        end
        
        % Get class time constants
        function [xE, dOfY, dt, dtE, xi, xf] = time_constants(obj)
            D = datevec(obj.date(1));
            dOfY = obj.date - datenum([D(1), 1, 1]);
            dt = obj.date(2:end) - obj.date(1:end-1);
            
            EplusTime = datenum(datevec(obj.mtr.Timestamp));
            D = datevec(EplusTime(1));
            EplusTime = EplusTime - datenum([D(1), 1,1]);
            
            xi = max(1,nnz((EplusTime)<=dOfY(1)));
            xf = nnz((EplusTime)<=dOfY(end));
            xE = EplusTime(xi:xf);
            if xi == 1
                dtE = xE - ([0;EplusTime(1:xf-1)]);
            else
                dtE = xE - (EplusTime(xi-1:xf-1));
            end
        end
        
        % Add plot elements
        function add_plot_elems(~, fig, figTitle)
            set(fig, 'name', figTitle)
            xlabel('Day of Year')
            ylabel('Electric Load (kW)')
            title(figTitle)
        end
        
        % Step plot
        function step_plot(~, varargin)
            xE = varargin{1};
            yE = varargin{2};
            if length(varargin) > 2
                lineSpec = {varargin{3:end}};
            else
                lineSpec = {};
            end
            stairs(xE, yE, lineSpec{:})
        end
        
        % Order open figures by figure number
        function order_figs(~)
            openFigs = findobj('Type', 'Figure');
            openFigNums = [];
            if ~isempty(openFigs)
                openFigNums = sort([openFigs(:).Number], 'descend');
            end
            for i = openFigNums
                figure(i)
            end
        end
    end
end
