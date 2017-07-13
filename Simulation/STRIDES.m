%% STRIDES %%
% This script gives the option to initialize or load a saved non-linear systems model
% Next it offers to linearize that model,or load a saved linearization. 
% Finally it is able to run a transient on either/both the linear or non-linear model
global Model_dir SimSettings WaitBar modelParam LinMod IterCount TagInf TagFinal Tags Inlet Outlet
%%clear variables
LinMod =[];
modelParam = [];
TagInf = [];
TagFinal = [];
Tags = [];
SimSettings = [];
Inlet = [];
Outlet = [];
%% Load a Plant
%%%%-- either load a plant and initialize it, or load a saved Plant & model Parameters
%%%%-- Working options are: SOFCstack, SOECstack, SOFCsystem, GasTurbine
J = center_menu('Non-linear model options','Initialize from system description','Load previously initialized plant','Skip to pre-loaded linear model');
if J ==1
    %the goal is to have a GUI replace these m-files
    ModelFiles=dir(fullfile(Model_dir,'Model Library','*.m'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(strrep(ModelFiles(i).name,'.m',''));
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        Plant = feval(modelName);
        WaitBar.Show = 1;
        BuildModel(Plant); %%Build & Initialize model
        J2 = center_menu('Save Model?','Yes','No');
        if J2 ==1
%             modelParam.SimSettings = SimSettings;
%             modelParam.Outlet = Outlet;
            [f,p]=uiputfile(fullfile(Model_dir,'Model Library','Saved Models',strcat(modelName,'.mat')),'Save Model As...');
            save([p f],'modelParam','SimSettings','Outlet')
        end
    end
elseif J ==2
    ModelFiles=dir(fullfile(Model_dir,'Model Library','Saved Models','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        modelName = list{s};
        modelName = strrep(modelName,'.mat','');
        load(fullfile(Model_dir,'Model Library','Saved Models',modelName));
%         SimSettings = modelParam.SimSettings;
%         Outlet = modelParam.Outlet;
    end
end
%% Create or load a set of linear models
nom = {};
if ~isempty(SimSettings)
    F = fieldnames(SimSettings);
    for i = 1:1:length(F)
        if strfind(F{i},'Nominal')
            nom(1,end+1) = F(i);
        end
    end
end
if isempty(nom) || J ==3
    J2 = 2; %load linearization
else J2 = center_menu('Linear model option','Create new linearization','Load previously linearized plant','Skip linearization');
end
if J2 ==1
    if length(nom)>1
        A = center_menu('Choose Nominal Condition to Linearize',nom);
        linParam = nom{A};
    else linParam = nom{1};
    end
    HSVtol = 0; %what should this be set to?
    Prompt = {'Minimum','Maximum','# of Linear Models'};
    DefaultVal = {num2str(0.2*SimSettings.(linParam)),num2str(SimSettings.(linParam)),'5'};
    A= inputdlg(Prompt,'Specify range and resolution of lineraization',1,DefaultVal);
    SetPoints = linspace(str2double(A(2)),str2double(A(1)),str2double(A(3)));
    CreateLinModel(SetPoints,HSVtol);
    J3 = center_menu('Save Linearized Model?','Yes','No');
    if J3 ==1
%         LinMod.SimSettings = SimSettings;
        [f,p]=uiputfile(fullfile(Model_dir,'Model Library','Saved Linearizations',strcat(modelName,'.mat')),'Save Linearized Model As...');
        save([p f],'LinMod','SimSettings','Outlet')
    end
elseif J2 ==2
    ModelFiles=dir(fullfile(Model_dir,'Model Library','Saved Linearizations','*.mat'));
    list = {};
    for i = 1:1:length(ModelFiles)
        list(i) = cellstr(ModelFiles(i).name);
    end
    [s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
    if OK~=1
        disp('Invalid selection. Exiting...')
    else
        load(fullfile(Model_dir,'Model Library','Saved Linearizations',list{s}));
%         if isempty(SimSettings)
%             SimSettings = LinMod.SimSettings;
%         end
%         if isfield(LinMod,'NominalPower')
%             SimSettings.NominalPower = LinMod.NominalPower;
%         end
    end
elseif J2 ==3
    LinMod =[];
end

if J2 ==3
    if J ==1 || J ==2
        J3 = center_menu('Simulation Options','Simulate non-linear response','Exit');
        if J3 == 2;
            J3 = 4;%exit
        end
    else J3 = 4; %no linear or non-linear model loaded, exit
    end
elseif J == 3
    J3 = 2; %simulate linear response
else
    J3 = center_menu('Simulation Options','Simulate non-linear response','Simulate linear response','Simulate both linear and non-linear response','Neither');
end
if J3 ~=4 
    %set up the transient
    %first identify any controller input (lookup functions), let user pick ones with a schedule
    %then get the variables from that function and allow the user to edit them
    A = (inputdlg('Test Duration (s)','Specify length of the transient simulation',1,{num2str(24*3600)}));
    SimSettings.RunTime = eval(A{1});
    if ~isempty(modelParam)
        controls = fieldnames(modelParam.Controls);
    elseif ~isempty(LinMod)
        controls = fieldnames(LinMod.Controls);
    end
    for i = 1:1:length(controls)
        if ~isempty(modelParam)
            Cont = modelParam.Controls.(controls{i});
        elseif ~isempty(LinMod)
            Cont = LinMod.Controls.(controls{i});
        end
        for j = 1:1:length(Cont.connections)
            r = strfind(Cont.connections{j},'.');
            if ~isempty(Cont.connections{j}) && ~isnumeric(Cont.connections{j}) && isempty(r) %must be a lookup function
                [globvar,Prompt,DefaultVal] = feval(Cont.connections{j},0,'loadparam'); 
                A = inputdlg(Prompt,strcat('Specify the transient imput parameters for the function',Cont.connections{j},'Any string will be evaluated, but must create vectors of equal length.'),1,DefaultVal);
                for k = 1:1:length(globvar)
                    SimSettings.(globvar{k}) = eval(A{k});
                end
            end
        end
        %% modify control terms
        nC = length(Cont.Gain);
        Prompt = {};
        DefaultVal = {};
        for j = 1:1:nC
            Prompt(end+1) = cellstr(strcat(Cont.PIdescription{j},'---Proportional'));
            Prompt(end+1) = cellstr(strcat(Cont.PIdescription{j},'---Integral'));
            DefaultVal(end+1) = cellstr(num2str(Cont.PropGain(j)));
            DefaultVal(end+1) = cellstr(num2str(Cont.Gain(j)));
        end
        A= inputdlg(Prompt,'Specify the controller parameters',1,DefaultVal);
        for j = 1:1:nC
            Cont.PropGain(j) = str2double(A(2*j-1));
            Cont.Gain(j) = str2double(A(2*j));
        end
        if ~isempty(modelParam)
            modelParam.Controls.(controls{i}) = Cont;
        end
        if ~isempty(LinMod)
            LinMod.Controls.(controls{i}) = Cont;
        end
        % %SOFCstack
        % Cont.Gain = [3e-3;1e-3;1e-2];
        % Cont.PropGain = [1;1;1];

        %SOFCsystem
        % Cont.Gain = [1e-2;1e-4;1e-2];
        % Cont.PropGain = [.5;.1;1];

        % %SOECstack
        % Cont.Gain = [3e-3;1e-3;1e-2];
        % Cont.PropGain = [1;1;1];

        % %GasTurbine
        % Cont.IntGain = [4e-4; 1e-2; 4e-2;];
        % Cont.PropGain = [8e-3; 5e-0; .75;];
    end 
end

%% Run a transient on non-linear model
if J3 ==1 || J3 == 3
    WaitBar.Show = 1;
    IterCount = 1; TagInf =[]; TagFinal =[];  WaitBar.Text = 'Running non-linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime], modelParam.IC); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    PlotSimulation(T,Y,1,0,1)% Plot the tags and scopes in Plant.Plot. The first option (after Y) is to plot any fuel cell or electrolyzer temperature profiles, the second option is to mave a video of the transient, the third is to plot compressor turbine and blower maps
end

%% Run transient on set of linear models
if J3 ==2 || J3 ==3
    % need to find better initial condition (won't always start at nominal power)
    IC =  [LinMod.Model{1}.X0;LinMod.Model{1}.UX0];
%     Tags.O = LinMod.Model{1}.Out0;
    Tags.U = LinMod.Model{1}.U0;
%     Tags.dX = zeros(length(LinMod.Model{1}.X0),1);
%     Tags.dUX = zeros(length(LinMod.Model{1}.UX0),1);
    IterCount = 1; TagInf =[]; TagFinal =[]; WaitBar.Show = 1; WaitBar.Text = 'Running linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    tic; [T, Y] = ode15s(@RunLinSystem, [0, SimSettings.RunTime], IC); disp(strcat('Time to run model:',num2str(toc),' seconds'));close(WaitBar.Handle);
    if J3 == 3
        PlotSimulation(T,Y,0,0,0) % have already plotted the maps, just adding non-linear response to linear response
    else
        PlotSimulation(T,Y,1,0,1)
    end
end