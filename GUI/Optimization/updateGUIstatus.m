function updateGUIstatus(handles,Solution,Date,Si)
global Plant Model_dir GENINDEX DateSim
if ~isempty(handles)
    if strcmp(get(handles.uipanelMain1,'Visible'),'on')
        nG = length(Plant.Generator);
        GenDisp = Solution.Dispatch(1:2,:);
        %% update status lights
        for i = 1:1:nG
            x = [];
            if Plant.Generator(i).Enabled
                if GenDisp(2,i)>0 && GenDisp(1,i)==0  %just turned on
                    [x,~] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
                end
                if (GenDisp(2,i)==0 && GenDisp(1,i)>0) || (Si==1 && GenDisp(2,i)==0) %just turned off or 1st time running
                    [x,~] = imread(fullfile(Model_dir,'GUI','Graphics','yellow.png'));
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
            if ~isempty(strfind(Plant.Generator(GENINDEX).Type,'Storage'))
                if strcmp(Plant.Generator(GENINDEX).Type,'Hydro Storage')
                    power = (GenDisp(2,GENINDEX));
                else
                    power = (GenDisp(1,GENINDEX)- GenDisp(2,GENINDEX))/Plant.optimoptions.Resolution*Plant.Generator(GENINDEX).QPform.Stor.DischEff;
                end
                set(handles.GenStatus1,'string', num2str(power));
                set(handles.GenStatus2text,'string','State-Of-Charge (%)');
                SOC = GenDisp(2,GENINDEX)/Plant.Generator(GENINDEX).QPform.Stor.UsableSize*100;
                set(handles.GenStatus2,'string', num2str(SOC));
            else
                set(handles.GenStatus1,'string', num2str(GenDisp(2,GENINDEX)));
                set(handles.GenStatus2text,'string','Efficiency (%)');
                skip = false;
                if ~isempty(Plant.Generator(GENINDEX).Output)
                    cap = Plant.Generator(GENINDEX).Output.Capacity*Plant.Generator(GENINDEX).Size;
                end
                if strcmp(Plant.Generator(GENINDEX).Type,'Electric Generator') || strcmp(Plant.Generator(GENINDEX).Type,'CHP Generator')
                    eff = Plant.Generator(GENINDEX).Output.Electricity;
                elseif strcmp(Plant.Generator(GENINDEX).Type,'Chiller') 
                    eff = Plant.Generator(GENINDEX).Output.Cooling;
                elseif strcmp(Plant.Generator(GENINDEX).Type,'Heater')
                    eff = Plant.Generator(GENINDEX).Output.Heat;    
                elseif strcmp(Plant.Generator(GENINDEX).Type,'Hydro')
                    skip = true;
                else skip = true;
                end
                if ~skip
                    Efficiency = interp1(cap,eff,GenDisp(2,GENINDEX))*100;
                    set(handles.GenStatus2,'string', num2str(Efficiency));
                end
            end
        end
        backSteps = min(Si-1,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
        History.Dispatch = Plant.Dispatch.GeneratorState(Si-backSteps:Si,:);
        History.LineFlows = Plant.Dispatch.LineFlows(Si-backSteps:Si,:);
        History.Buildings = Plant.Dispatch.Buildings(Si-backSteps:Si,:);
        History.hydroSOC = Plant.Dispatch.hydroSOC(Si-backSteps:Si,:);
        HistoryTime = Plant.Dispatch.Timestamp(Si-backSteps:Si,:);
        plotDispatch(handles,Date,Solution,HistoryTime,History)
        set(handles.sliderStartDate,'Value',DateSim - Plant.Dispatch.Timestamp(1))
    elseif strcmp(get(handles.uipanelMain2,'Visible'),'on')
        Plant.Market.MarginCost = MarginalCapacityCost(Solution.Dispatch,Date);
        plotMarginalCapacityCost(handles)
    end
end
end%ends function updateGUIstatus