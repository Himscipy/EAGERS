function update_gui_status(handles,solution,date,s_i,gen,options,dispatch)
global GENINDEX
dir = strrep(which('update_gui_status.m'),fullfile('GUI','Optimization','update_gui_status.m'),'');
if ~isempty(handles)
    if strcmp(get(handles.uipanelMain1,'Visible'),'on')
        nG = length(gen);
        GenDisp = solution.Dispatch(1:2,:);
        %% update status lights
        for i = 1:1:nG
            x = [];
            if gen(i).Enabled
                if GenDisp(2,i)>0 && GenDisp(1,i)==0  %just turned on
                    [x,~] = imread(fullfile(dir,'GUI','Graphics','green.png'));
                end
                if (GenDisp(2,i)==0 && GenDisp(1,i)>0) || (s_i==1 && GenDisp(2,i)==0) %just turned off or 1st time running
                    [x,~] = imread(fullfile(dir,'GUI','Graphics','yellow.png'));
                end
            end
            if ~isempty(x)
                if license('test','Image Processing Toolbox')
                    set(handles.Switch,'Units','pixels');
                    pos1 = get(handles.Switch,'Position');
                    set(handles.Switch,'Units','characters');
                    pos2 = get(handles.Switch,'Position');
                    x = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
                end
                set(handles.(strcat('GeneratorStat_',num2str(i))),'cdata', x);
            end
        end

        %% Update status of selected generator in lower left box
        if ~isempty(GENINDEX)
            gen_i = gen(GENINDEX);
            if ~isempty(strfind(gen_i.Type,'Storage'))
                if strcmp(gen_i.Type,'Hydro Storage')
                    power = (GenDisp(2,GENINDEX));
                else
                    power = (GenDisp(1,GENINDEX)- GenDisp(2,GENINDEX))/options.Resolution*gen_i.QPform.Stor.DischEff;
                end
                set(handles.GenStatus1,'string', num2str(power));
                set(handles.GenStatus2text,'string','State-Of-Charge (%)');
                SOC = GenDisp(2,GENINDEX)/gen_i.QPform.Stor.UsableSize*100;
                set(handles.GenStatus2,'string', num2str(SOC));
            else
                set(handles.GenStatus1,'string', num2str(GenDisp(2,GENINDEX)));
                set(handles.GenStatus2text,'string','Efficiency (%)');
                skip = false;
                if ~isempty(gen_i.Output)
                    cap = gen_i.Output.Capacity*gen_i.Size;
                end
                if strcmp(gen_i.Type,'Electric Generator') || strcmp(gen_i.Type,'CHP Generator')
                    eff = gen_i.Output.Electricity;
                elseif strcmp(gen_i.Type,'Chiller') 
                    eff = gen_i.Output.Cooling;
                elseif strcmp(gen_i.Type,'Heater')
                    eff = gen_i.Output.Heat;    
                elseif strcmp(gen_i.Type,'Hydro')
                    skip = true;
                else
                    skip = true;
                end
                if ~skip
                    Efficiency = interp1(cap,eff,GenDisp(2,GENINDEX))*100;
                    set(handles.GenStatus2,'string', num2str(Efficiency));
                end
            end
        end
        backSteps = min(s_i-1,options.Horizon/options.Resolution);
        History.Dispatch = dispatch.GeneratorState(s_i-backSteps:s_i,:);
        History.LineFlows = dispatch.LineFlows(s_i-backSteps:s_i,:);
        History.Buildings = dispatch.Buildings(s_i-backSteps:s_i,:);
        History.hydroSOC = dispatch.hydroSOC(s_i-backSteps:s_i,:);
        HistoryTime = dispatch.Timestamp(s_i-backSteps:s_i,:);
        plotDispatch(handles,date,solution,HistoryTime,History)
        set(handles.sliderStartDate,'Value',date(1) - dispatch.Timestamp(1))
    elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
        market_margin_cost = marginal_cost(gen,solution.Dispatch,date);
        plotMarginalCapacityCost(handles,market_margin_cost,date)
    end
end
end%ends function updateGUIstatus