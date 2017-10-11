function dX = RunLinSystem(t,X)
%%RunLinSystem: Simulates the transinet response by interpolating a set 
%%of linear time invarient models that have the same states
%%The model is linear, but it can be connected to any controller function
%%The inputs are time (in seconds) and X, the vector of states
%%This is solving a model of the form:
%%dX = A*x+B*u
%%Y = C*x+D*u
%%with a controller u = fcn(Y,inputs or forcing functions)
global LinMod IterCount TagInf Tags WaitBar SimSettings Jcount

%%this section sets up the counters
nS = length(LinMod.Model{1}.A);%number of model states
if isempty(TagInf) %create TagInf Fields
    TagInf.Time(1) =0;
    for i = 1:1:length(LinMod.Scope)
        r = strfind(LinMod.Scope{i},'.');
        TagInf.(LinMod.Scope{i}(1:r-1)).(LinMod.Scope{i}(r+1:end)) = Tags.(LinMod.Scope{i}(1:r-1)).(LinMod.Scope{i}(r+1:end));
    end
    Jcount = 0;
end
if t==0 || t == TagInf.Time(IterCount)
    Jcount = Jcount+1;
elseif t>TagInf.Time(IterCount)%+1e-5
    IterCount =IterCount+1;
    Jcount = nS;
elseif IterCount>1 && t<TagInf.Time(IterCount)
    Jcount = 1;
    while t<TagInf.Time(IterCount-1)
        IterCount = IterCount-1;
    end
end
TagInf.Time(IterCount,1) = t;

%%This section sets up the interpolation between the models
ControlInlet = ModelOutlet(t,LinMod.Model{1}.Out0);%will update inlets, specifically the controller target that we linearized around
Target = LinMod.LinearizedTarget;
r = strfind(Target,'.');
InterpVal = ControlInlet.(Target(1:r-1)).(Target(r+1:end)); 
c= [0 0 0 0];
n = length(LinMod.InterpVec);
if LinMod.InterpVec(1) > LinMod.InterpVec(2) %values of InterpVec are decreasing
    k = nnz(LinMod.InterpVec > InterpVal);
else
    k = nnz(LinMod.InterpVec < InterpVal);
end
if k == n
    nModel(4) = n;
    c(4) = 1; 
elseif k==0
    nModel(1) = 1;
    c(1) = 1; 
elseif k == n-1
    nModel(4) = k+1;
    nModel(3) = k;
    c(4) = (InterpVal - LinMod.InterpVec(n-1))/(LinMod.InterpVec(n) - LinMod.InterpVec(n-1));
    c(3) = (InterpVal - LinMod.InterpVec(n))/(LinMod.InterpVec(n-1) - LinMod.InterpVec(n));
elseif k == 1
    nModel(2) = k+1;
    nModel(1) = k;
    c(2) = (InterpVal - LinMod.InterpVec(1))/(LinMod.InterpVec(2) - LinMod.InterpVec(1));
    c(1) = (InterpVal - LinMod.InterpVec(2))/(LinMod.InterpVec(1) - LinMod.InterpVec(2));
else 
    nModel(4) = k+2;
    nModel(3) = k+1;
    nModel(2) = k;
    nModel(1) = k-1;
    Vec = [LinMod.InterpVec(k-1),LinMod.InterpVec(k),LinMod.InterpVec(k+1),LinMod.InterpVec(k+2)];
    c = ones(1,4);
    for i = 1:4
        for j = 1:4
            if i ~= j
                c(i) = c(i)*(InterpVal - Vec(j))/(Vec(i) - Vec(j));
            end
        end
    end
end

%% linear model + controller to find dX and dUX
ControlStates = X(nS+1:end,1); %seperate the controller states (at end of vector)
X = X(1:nS,1); %reduce to just the model states

U = LinMod.Input; %controller outputs from last iteration
for j = 1:1:3
    dX = 0;
    Y = 0;
    for i = 1:1:4
        if c(i)>0
            dX = dX + c(i)*(LinMod.Model{nModel(i)}.A*(X-LinMod.Model{nModel(i)}.X0) + LinMod.Model{nModel(i)}.B*(U-LinMod.Model{nModel(i)}.U0));%calculate dX from X's and U's 
            Y = Y + c(i)*(LinMod.Model{nModel(i)}.Out0 + (LinMod.Model{nModel(i)}.C*(X-LinMod.Model{nModel(i)}.X0) + LinMod.Model{nModel(i)}.D*(U-LinMod.Model{nModel(i)}.U0)));%calculate O from X's and U's
        end
    end
    [U,dU] = Controller(t,Y,ControlStates);
end

LinMod.Input = U;
dX(end+1:end+length(ControlStates)) = dU;

%%put infor into dialog box to show status
if t>0 && WaitBar.Show == 1 && Jcount==length(X) && isfield(LinMod,'Scope')
    n = length(LinMod.Scope);
    dt = TagInf.Time(IterCount,1) - TagInf.Time(IterCount-1,1);
    points = nnz(TagInf.Time>(TagInf.Time(IterCount,1)-100*dt));
    for i = 1:1:n
        tagName = LinMod.Scope{i};
        r = strfind(tagName,'.');
        block = tagName(1:r-1);
        tag = tagName(r+1:end);
        if isfield(TagInf,block) && isfield(TagInf.(block),tag)
            figure(i)
            plot(TagInf.Time(IterCount-points+1:IterCount,1),TagInf.(block).(tag)(IterCount-points+1:IterCount,:))
            ylabel(tagName);
            xlabel('Time in seconds');
        end
    end
end
%% update waitbar
if t ==0
    Text = 'calculating Jacobian';
    x = Jcount/length(X);
elseif IterCount>1 && Jcount<length(X)
    Text = 're-calculating Jacobian';
    x = Jcount/length(X);
else
    if strcmp(WaitBar.Text,'Running linear model with transient')
        x = t/SimSettings.RunTime;
    else x = max(0,(log(t/1e-4))/(log(SimSettings.RunTime/1e-4)));
    end
    Text = {WaitBar.Text;strcat('    Time = ', num2str(t))};
end
if WaitBar.Show == 1
    waitbar(x,WaitBar.Handle,Text);
end
end%Ends function RunLinSystem

function [U,dU] = Controller(t,ModelOut,ControlStates)
% what to do if the controller outputs and # of states don't align?
global LinMod Tags TagInf IterCount
ControlInlet = ModelOutlet(t,ModelOut);
list = fieldnames(LinMod.Controls);
n = 0;
for i = 1:1:length(list)
    block = list{i};
    s = length(LinMod.Controls.(block).Scale);
    Yi = ControlStates(n+1:n+s,1).*LinMod.Controls.(block).Scale;
    ControllerOut.(block) = feval(LinMod.Controls.(block).type,t,Yi,ControlInlet.(block),LinMod.Controls.(block),'Outlet');
    dU(n+1:n+s,1) = feval(LinMod.Controls.(block).type,t,Yi,ControlInlet.(block),LinMod.Controls.(block),'dY');
    dU(n+1:n+s,1) = dU(n+1:n+s,1)./LinMod.Controls.(block).Scale;
    n = n+s;
end
U = ModelInlet(ControllerOut);  
if isfield(LinMod.Controls.(block),'TagInf')
    tagNames = LinMod.Controls.(block).TagInf;
    for i = 1:1:length(tagNames)
        if isnumeric(Tags.(block).(tagNames{i}))
            TagInf.(block).(tagNames{i})(IterCount,:)=Tags.(block).(tagNames{i});
        else
            f = fieldnames(Tags.(block).(tagNames{i}));
            for j = 1:1:length(f)
                TagInf.(block).(tagNames{i}).(f{j})(IterCount,:)=Tags.(block).(tagNames{i}).(f{j}); 
            end
        end
    end
end
end%Ends function Controller


function U0 = ModelInlet(ControllerOut)
global LinMod Tags
%this creates a vector (U0) of values from the controller and any lookup functions or tags that feed into the linear model
components = fieldnames(LinMod.Components);
controls = fieldnames(LinMod.Controls);
U0 = []; 
for k = 1:1:length(components)
    list = LinMod.Components.(components{k}).InletPorts; %only field should be ports with a connection to a block, lookup function or tag
    for i = 1:1:length(list)
        port = list{i};
        A = [];
        BlockPort = char(LinMod.Components.(components{k}).(port).connected);
        r = strfind(BlockPort,'.');
        if ~isempty(r)
            connectedBlock = BlockPort(1:r-1);
            if strcmp(connectedBlock,'Tags')
                connectedBlock = BlockPort(r(1)+1:r(2)-1);
                connectedPort = BlockPort(r(2)+1:end);
                A = Tags.(connectedBlock).(connectedPort);
            elseif ismember(connectedBlock,controls)
                connectedPort = BlockPort(r+1:end);
                A = ControllerOut.(connectedBlock).(connectedPort);
            end
        elseif ~isempty(BlockPort)
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

function ControlInlet = ModelOutlet(t,Out)
global LinMod TagInf IterCount
%this creates a vector (Out0) of values for the controller inlet and any scopes, these will be the outputs of the linear model
controls = fieldnames(LinMod.Controls);
components = fieldnames(LinMod.Components);
n = 0;
for k = 1:1:length(controls)
    block = controls{k};
    list = LinMod.Controls.(block).InletPorts; 
    for i = 1:1:length(list)
        port = list{i};
        A = []; %In this function A is a placeholder identifying how many Out values get sent to that controller input port
        if ~isempty(LinMod.Controls.(block).(port).connected)%inlet connected to another outlet
            BlockPort = char(LinMod.Controls.(block).(port).connected);
            if ~isnumeric(BlockPort) %ignore constant inlet values
                r = strfind(BlockPort,'.');
                if ~isempty(r)
                    connectedBlock = BlockPort(1:r-1);
                    if strcmp(connectedBlock,'Tags') || ismember(connectedBlock,components)
                        A = LinMod.Controls.(block).(port).IC; %port is connected to the model or a tag
                    end
                else
                    ControlInlet.(block).(port) = feval(BlockPort,t);%finds value of look-up function at time = 0
                end
            end
        else
            ControlInlet.(block).(port) = LinMod.Controls.(block).(port).IC;%not connected use constant value
        end
        if ~isempty(A)
            if isstruct(A)
                F = fieldnames(A);
                for j = 1:1:length(F)
                    s = length(A.(F{j}));
                    ControlInlet.(block).(port).(F{j}) = Out(n+1:n+s,1);
                    n = n+s;
                end
            else
                s = length(A);
                ControlInlet.(block).(port) = Out(n+1:n+s,1);
                n = n+s;
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
    A = TagInf.(connectedBlock).(connectedTag)(1,:);
    if isstruct(A)
        F = fieldnames(A);
        for j = 1:1:length(F)
            s = length(A.(F{j}));
            TagInf.(connectedBlock).(connectedTag).(F{j})(IterCount,:) = Out(n+1:n+s,1)';
            n = n+s;
        end
    else
        s = length(A);
        TagInf.(connectedBlock).(connectedTag)(IterCount,:) = Out(n+1:n+s,1)';
        n = n+s;
    end
end
end%Ends function ModelOutlet