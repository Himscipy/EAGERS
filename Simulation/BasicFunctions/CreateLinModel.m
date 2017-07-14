function CreateLinModel
%builds linear time independent model around a steady state point for a model
global  modelParam WaitBar SimSettings IterCount TagInf TagFinal LinMod 
controls = fieldnames(modelParam.Controls);
for i = 1:length(controls)
    LinMod.Controls.(controls{i}) = modelParam.Controls.(controls{i});
end
ControlStates = [];%find controller states
for i = 1:length(controls)
    ControlStates((end+1):(end+length(modelParam.Controls.(controls{i}).States))) = modelParam.Controls.(controls{i}).States;
end
components = fieldnames(modelParam.Components);
for i = 1:length(components)
    list2 = modelParam.Components.(components{i}).InletPorts;
    for j = 1:1:length(list2)
        port = list2{j};
        if ~isempty(modelParam.Components.(components{i}).(port).connected)
            LinMod.Components.(components{i}).(port) = modelParam.Components.(components{i}).(port).connected;%connected to another block
%         else LinMod.Components.(components{i}).(port) = modelParam.Components.(components{i}).(port).IC;%constant value, won't be part of linear model input since it is fixed
        end
    end
end
if isfield(modelParam,'Scope')
    LinMod.Scope = modelParam.Scope;
end
LinMod.UpperBound = modelParam.UpperBound;
LinMod.LowerBound = modelParam.LowerBound;
LinMod.IC = modelParam.IC;
LinMod.Input = ModelInlet; %initial model input
Y0 =  LinMod.IC;
controls = fieldnames(LinMod.Controls);
i = 0;
while i<length(controls)
    i = i+1;
    targetList = LinMod.Controls.(controls{i}).TargetDescription;
    targetList(end+1) = {'none of these'};
    L1 = center_menu('Which controller target would you like to linearize around?',targetList);
    if L1~=length(targetList)
        LinMod.LinearizedTarget = strcat(controls{i},'.Target',num2str(L1));
        targetSource = LinMod.Controls.(controls{i}).connections{L1};
        [globvar,Prompt,DefaultVal] = feval(targetSource,0,'loadparam'); 
        SimSettings.RunTime = str2double(inputdlg('Time to ensure steady-state (hr)','Simulation duration',1,{'4'}))*3600;
        L2 = center_menu('Select linearization variable',Prompt);
        Prompt = {'Minimum','Maximum','# of Linear Models'};
        DefaultVal = {num2str(0.2*SimSettings.(globvar{L2})(1)),num2str(SimSettings.(globvar{L2})(1)),'5'};
        A = inputdlg(Prompt,'Specify range and resolution of lineraization',1,DefaultVal);
        SetPoints = linspace(str2double(A(2)),str2double(A(1)),str2double(A(3)));
        LinMod.InterpVec = SetPoints;
        i = length(controls);
    end
end
%make sure all time variables in SimSettings reach to RunTime (note these are in hours)
F = fieldnames(SimSettings);
for i = 1:1:length(F)
    if strcmp(F{i},'RunTime')
        %do nothing
    elseif ~isempty(strfind(F{i},'Time')) && length(SimSettings.(F{i}))>1
        SimSettings.(F{i}) = [0, SimSettings.RunTime/3600];
    elseif length(SimSettings.(F{i}))>1
        SimSettings.(F{i}) = [SimSettings.(F{i})(1), SimSettings.(F{i})(1)];%hold inputs constant
    end
end

for n = 1:length(SetPoints)
    SimSettings.(globvar{L2}) = [SimSettings.(globvar{L2})(end),SetPoints(n)];
    IterCount = 1; TagInf =[]; TagFinal =[]; LinMod.Interupt = []; WaitBar.Show = 0;
%     WaitBar.Show = 1;  WaitBar.Text = 'Running non-linear model with transient';WaitBar.Handle =waitbar(0,WaitBar.Text);
    [~, Y] = ode15s(@RunBlocks, [0, SimSettings.RunTime],Y0);
%     close(WaitBar.Handle);
%     WaitBar.Show = 0;
    Y0 = Y(end,:)';%states at steady state
    U0 = ModelInlet;% convert model inputs to vector (if coming from controller, tag, or lookup function)
    Out0 = ModelOutlet;% convert model outputs to vector (if going to controller, or scope)

    %% find effect of state perturbations in the model portion
    ModelStates = (1:1:length(Y0))';
    ModelStates(ControlStates) = [];
    Y1 = Y0(ModelStates);
    
    A = zeros(length(Y1));
    C = zeros(length(Out0));
    Perturbation = eye(length(Y0))*1e-9;
    for i = 1:length(Y1)%for all non-ignored states
        epsY = 2e-9;
        if (Y1(i)<0 && Y1(i)>-1e-9) || Y1(i)==LinMod.UpperBound(ModelStates(i)) % avoid switching signs on states or passing an absolute threshold
            dYplus = 0*Y1;
            OutPlus = Out0;
            epsY = 1e-9;
        else
            dYplus = RunBlocks(SimSettings.RunTime,Y0+Perturbation(:,ModelStates(i)));%record resulting dY's
            dYplus(ControlStates)=[];%ignore dY's from ignored states
            OutPlus = ModelOutlet;
        end
        if (Y1(i)>=0 && Y1(i)<1e-9) || Y1(i)==LinMod.LowerBound(ModelStates(i)) % avoid switching signs on states or passing an absolute threshold
            dYminus = 0*Y1;
            OutMinus = Out0;
            epsY = 1e-9;
        else
            dYminus = RunBlocks(SimSettings.RunTime,Y0-Perturbation(:,ModelStates(i)));
            dYminus(ControlStates) = [];
            OutMinus = ModelOutlet;
        end
        A(:,i) = (dYplus-dYminus)/epsY;
        C(:,i) = (OutPlus-OutMinus)/epsY;
    end

    %% find effect of input perturbations on dY and Outputs
    B = zeros(length(Y1),length(U0));
    D = zeros(length(Out0),length(U0));
    Perturbation = eye(length(U0))*1e-9;
    for i = 1:length(U0)%for all inputs
        epsU = 2e-9;
        if (U0(i)<0 && U0(i)>-1e-9) || U0(i)==LinMod.UpperBound(ControlStates(i)) % avoid switching signs on states or passing an absolute threshold
            dYplus = 0*Y1;
            OutPlus = Out0;
            epsU = 1e-9;
        else
            OverideInput(U0 + Perturbation(:,i),i);
            dYplus = RunBlocks(SimSettings.RunTime,Y0);%record resulting dY's
            dYplus(ControlStates)=[];%ignore dY's from ignored states
            OutPlus = ModelOutlet;
        end
        if (U0(i)>0 && U0(i)<1e-9) || U0(i)==LinMod.LowerBound(ControlStates(i)) % avoid switching signs on states or passing an absolute threshold
            dYminus = 0*Y1;
            OutMinus = Out0;
            epsU = 1e-9;
        else
            OverideInput(U0 - Perturbation(:,i),i);
            dYminus = RunBlocks(SimSettings.RunTime,Y0);
            dYminus(ControlStates) = [];
            OutMinus = ModelOutlet;
        end

        B(:,i) = (dYplus - dYminus)/epsU;
        D(:,i) = (OutPlus-OutMinus)/epsU;
    end

    Model.A = A;
    Model.B = B;
    Model.C = C;
    Model.D = D;

    Model.U0 = U0;
    Model.Out0 = Out0;
    Model.UX0 = Y0(ControlStates);
    Model.X = Y1;
    LinMod.Options = [];
    LinMod.Model{n} = Model;
    disp(strcat('Created linear model_',num2str(n),'_of_',num2str(length(SetPoints))))
end

%create balanced realization and transform matrix from middle setpoint
n = round(length(SetPoints)/2);
Sys = ss(LinMod.Model{n}.A,LinMod.Model{n}.B,LinMod.Model{n}.C,LinMod.Model{n}.D);% convert to model
[Sys,g,T,Ti] = balreal(Sys);
LinMod.Model{n}.Sys = Sys;
LinMod.Model{n}.HSV = g;
LinMod.Model{n}.T = T;
LinMod.Model{n}.Ti = Ti;

%removes states using balred and hankel singular values
HSVtol = 0; %what should this be set to?
keep = (1:length(LinMod.Model{1}.A));
keep(LinMod.Model{n}.HSV < HSVtol) = [];
for i = 1:length(SetPoints)
    A = LinMod.Model{n}.T*LinMod.Model{i}.A*LinMod.Model{n}.Ti;
    B = LinMod.Model{n}.T*LinMod.Model{i}.B;
    C = LinMod.Model{i}.C*LinMod.Model{n}.Ti;
    X0 = LinMod.Model{n}.T*LinMod.Model{i}.X;

    LinMod.Model{i}.A = A(keep,keep);
    LinMod.Model{i}.B = B(keep,:);
    LinMod.Model{i}.C = C(:,keep);
    LinMod.Model{i}.X0 = X0(keep);
end
end%Ends function CreateLinModel

function OverideInput(Unew,interupt)
global LinMod Outlet Tags
%this creates a vector (U0) of values from the controller and any lookup functions or tags that feed into the linear model
components = fieldnames(LinMod.Components);
controls = fieldnames(LinMod.Controls);
n = 0;
k = 0;
while n<interupt
    k = k+1;
    block = components{k};
    list = fieldnames(LinMod.Components.(block)); %only field should be ports with a connection to a block, lookup function or tag
    i = 0;
    while i<length(list)
        i = i+1;
        port = list{i};
        A = [];
        BlockPort = char(LinMod.Components.(block).(port));
        r = strfind(BlockPort,'.');
        if ~isempty(r)
            connectedBlock = BlockPort(1:r-1);
            if strcmp(connectedBlock,'Tags')
                connectedBlock = BlockPort(r(1)+1:r(2)-1);
                connectedPort = BlockPort(r(2)+1:end);
                A = Tags.(connectedBlock).(connectedPort);
            elseif ismember(connectedBlock,controls)
                connectedPort = BlockPort(r+1:end);
                A = Outlet.(connectedBlock).(connectedPort);
            end
        else
            A = feval(BlockPort,0);%finds value of look-up function at time = 0
        end
        if ~isempty(A)
            if isstruct(A)
                F = fieldnames(A);
                j = 0;
                while j<length(F)
                    j = j+1;
                    s = length(A.(F{j}));
                    n = n+s;
                    if n>=interupt
                        LinMod.Interupt.block = block;
                        LinMod.Interupt.port = port;
                        LinMod.Interupt.struct = F{j};
                        LinMod.Interupt.index = interupt - (n-s);
                        LinMod.Interupt.value = Unew(interupt);
                        j = length(F);%exit while loop
                        i = length(list);%exit inner while loop
                    end
                end
            else
                s = length(A);
                n = n+s;
                if n>=interupt
                    LinMod.Interupt.block = block;
                    LinMod.Interupt.port = port;
                    LinMod.Interupt.struct = {};
                    LinMod.Interupt.index = interupt - (n-s);
                    LinMod.Interupt.value = Unew(interupt);
                    i = length(list);%exit inner while loop
                end
            end
        end
    end
end
end%Ends function OverRideInput

function U0 = ModelInlet
global LinMod Outlet Tags
%this creates a vector (U0) of values from the controller and any lookup functions or tags that feed into the linear model
components = fieldnames(LinMod.Components);
controls = fieldnames(LinMod.Controls);
U0 = []; 
for k = 1:1:length(components)
    list = fieldnames(LinMod.Components.(components{k})); %only field should be ports with a connection to a block, lookup function or tag
    for i = 1:1:length(list)
        port = list{i};
        A = [];
        BlockPort = char(LinMod.Components.(components{k}).(port));
        r = strfind(BlockPort,'.');
        if ~isempty(r)
            connectedBlock = BlockPort(1:r-1);
            if strcmp(connectedBlock,'Tags')
                connectedBlock = BlockPort(r(1)+1:r(2)-1);
                connectedPort = BlockPort(r(2)+1:end);
                A = Tags.(connectedBlock).(connectedPort);
            elseif ismember(connectedBlock,controls)
                connectedPort = BlockPort(r+1:end);
                A = Outlet.(connectedBlock).(connectedPort);
            end
        else
            A = feval(BlockPort,0);%finds value of look-up function at time = 0
        end
        if ~isempty(A)
            if isstruct(A)
                F = fieldnames(A);
                for j = 1:1:length(F)
                    s = length(A.(F{j}));
                    U0(end+1:end+s,1) = A.(F{j});
                end
            else
                s = length(A);
                U0(end+1:end+s,1) = A;
            end
        end
    end
end
end%Ends function ModelInlet

function Out0 = ModelOutlet
global LinMod Outlet Tags
%this creates a vector (Out0) of values for the controller inlet and any scopes, these will be the outputs of the linear model
controls = fieldnames(LinMod.Controls);
components = fieldnames(LinMod.Components);
Out0 = []; 
for k = 1:1:length(controls)
    list = LinMod.Controls.(controls{k}).InletPorts; 
    for i = 1:1:length(list)
        port = list{i};
        A = [];
        if ~isempty(LinMod.Controls.(controls{k}).(port).connected)%inlet connected to another outlet
            BlockPort = char(LinMod.Controls.(controls{k}).(port).connected);
            if ~isnumeric(BlockPort) %ignore constant inlet values
                r = strfind(BlockPort,'.');
                if ~isempty(r)
                    connectedBlock = BlockPort(1:r-1);
                    if strcmp(connectedBlock,'Tags')
                        connectedBlock = BlockPort(r(1)+1:r(2)-1);
                        connectedTag = BlockPort(r(2)+1:end);
                        A = Tags.(connectedBlock).(connectedTag);
                    elseif ismember(connectedBlock,components)
                        connectedPort = BlockPort(r+1:end);
                        A = Outlet.(connectedBlock).(connectedPort);
                    end
                end
            end
        end
        if ~isempty(A)
            if isstruct(A)
                F = fieldnames(A);
                for j = 1:1:length(F)
                    s = length(A.(F{j}));
                    Out0(end+1:end+s,1) = A.(F{j});
                end
            else
                s = length(A);
                Out0(end+1:end+s,1) = A;
            end
        end
    end
end
%append with any scopes
for i = 1:1:length(LinMod.Scope)
    BlockPort = LinMod.Scope{i};
    r = strfind(BlockPort,'.');
    connectedBlock = BlockPort(1:r-1);
    connectedTag = BlockPort(r+1:end);
    A = Tags.(connectedBlock).(connectedTag);
    if isstruct(A)
        F = fieldnames(A);
        for j = 1:1:length(F)
            s = length(A.(F{j}));
            Out0(end+1:end+s,1) = A.(F{j});
        end
    else
        s = length(A);
        Out0(end+1:end+s,1) = A;
    end
end
end%Ends function ModelOutlet