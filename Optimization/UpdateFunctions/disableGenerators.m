function QP = disableGenerators(QP,Locked,Enabled)
% QP are the updated QP matrices
% Organize is the record of which indices are associated with each generator
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Enabled is a variable to include each generator in the optimization or not.
[m,n] = size(QP.organize);
nS = m-1;
nG = length(QP.Organize.Dispatchable);
nB = length(QP.Organize.Building.r);
nL = n - nG - nB;
if ~isempty(Locked)
    Enabled = true(nG,1);
    for j = 1:1:nG
        if all(Locked(:,j)==0)
            Enabled(j) = false;
        end
    end
end

%add in constant electric loads from chillers/pumps in FitB
if isfield(QP,'constDemand') 
    Outs = fieldnames(QP.constDemand);
    for s = 1:1:length(Outs)
        if any(QP.constDemand.(Outs{s}).req>0)
            if ~isempty(Locked)
                for j = 1:1:nG
                    QP.constDemand.(Outs{s}).load(:,j) = QP.constDemand.(Outs{s}).load(:,j).*Locked(2:end,j);
                end
            elseif ~isempty(Enabled)
                QP.constDemand.(Outs{s}).load(:,~Enabled) = 0;
            end
            QP.beq(QP.constDemand.(Outs{s}).req) = QP.beq(QP.constDemand.(Outs{s}).req) + sum(QP.constDemand.(Outs{s}).load,2);
        end
    end
end
    
%remove disabled generators from optimization matrices
rmv = false;
for i = 1:1:nG
    if ~Enabled(i)  && ~isempty(QP.Organize.States{i})
        rmv = true;  %there is a disabled generator with states to remove
        break
    end
end

[r,~] = size(QP.A);
[req,xL] = size(QP.Aeq);
rkeep = linspace(1,r,r)';
reqkeep = linspace(1,req,req)';
xkeep = linspace(1,xL,xL)';
if rmv    
    %% remove entire generators
    if nS == 0 %single time step
        s_rm = 0;
        r_rm = 0;
        req_rm = 0;
        for i = 1:1:n
            if ~isempty(QP.Organize.States{i})
                if i<=nG && ~Enabled(i)%remove disabled generator
                    QP.organize{1,i} = [];
                    states = QP.Organize.States{i};
                    xkeep(states) = 0; %removes states associated with this generator.
                    QP.Organize.States(i) = {[]}; %empty
                    s_rm = s_rm+length(states);
                    if QP.Organize.SpinReserveStates(i)>0
                        xkeep(QP.Organize.SpinReserveStates(i)) = 0; %removes state associated with this spinning reserve for this generator.
                        QP.Organize.SpinReserveStates(i) = 0;
                        s_rm = s_rm+1;
                    end
                    
                    if QP.Organize.Equalities(i,1)>0
                        eq = QP.Organize.Equalities(i,:);
                        reqkeep(eq(1):eq(1)+eq(2)-1) = 0;%equality constraints
                        QP.Organize.Equalities(i,1:2) = [0,0];
                        req_rm = req_rm + eq(2);
                    end
                    
                    if QP.Organize.SpinRow(i)>0
                        reqkeep(QP.Organize.SpinRow(i):QP.Organize.SpinRow(i)+1) = 0;% spinning reserve equality constraint 1 & 2
                        req_rm = req_rm + 2;
                    end
                    if QP.Organize.Inequalities(i)>0
                        rkeep(QP.Organize.Inequalities(i)) = 0;%1 inequality constraint for storage charging state
                        r_rm = r_rm + 1;
                    end
                        
                else
                    %don't remove, but change state #'s
                    QP.Organize.States(i) = {QP.Organize.States{i}-s_rm}; %empty
                    if i<=length(QP.Organize.SpinReserveStates)
                        if QP.Organize.SpinReserveStates(i)>0
                            QP.Organize.SpinReserveStates(i) = QP.Organize.SpinReserveStates(i) - s_rm;
                        end
                    end
                    if ~isempty(QP.organize{1,i})
                        QP.organize{1,i} = QP.organize{1,i} -s_rm; %
                    end
                    if i<=nG
                        if QP.Organize.Equalities(i,1)>0
                            QP.Organize.Equalities(i,1) = QP.Organize.Equalities(i,1) - req_rm;
                        end
                        if QP.Organize.SpinRow(i)>0
                            QP.Organize.SpinRow(i) = QP.Organize.SpinRow(i) - req_rm;
                        end
                        if QP.Organize.Inequalities(i)>0
                            QP.Organize.Inequalities(i) = QP.Organize.Inequalities(i) - r_rm;
                        end
                    end
                end
            end
        end
        QP.Organize.t1States = QP.Organize.t1States - s_rm;
        QP.Organize.t1ineq = QP.Organize.t1ineq - r_rm;
        QP.Organize.t1Balances = QP.Organize.t1Balances - req_rm;
    else %multi-time step
        ic_rm = 0;%initial conditions removed
        if isfield(QP.Organize,'IC')
            for i = 1:1:n
                if QP.Organize.IC(i)>0
                    if i<=nG && ~Enabled(i)% || (~isempty(Locked) && Locked(1,i)==0))
                        xkeep(QP.Organize.IC(i)) = 0;
                        reqkeep(QP.Organize.IC(i)) = 0; %initial condition constraints
                        QP.Organize.IC(i) = 0;
                        ic_rm = ic_rm+1;
                    else
                        QP.Organize.IC(i) = QP.Organize.IC(i)-ic_rm;
                    end
                end
            end
        end
        s_rm = ic_rm;
        r_rm = 0;
        req_rm = ic_rm;
        for i = 1:1:n
            if ~isempty(QP.Organize.States{i})
                if i<=nG && ~Enabled(i)%remove disabled generator
                    QP.organize{2,i} = [];

                    states = QP.Organize.States{i};
                    for j = 1:1:length(states)
                        xkeep(states(j):QP.Organize.t1States:(nS-1)*QP.Organize.t1States+states(j)) = 0; %removes states associated with this generator.
                    end
                    QP.Organize.States(i) = {[]}; %empty
                    s_rm = s_rm+length(states);

                    if QP.Organize.SpinReserveStates(i)>0
                        xkeep(QP.Organize.SpinReserveStates(i):QP.Organize.t1States:(nS-1)*QP.Organize.t1States+QP.Organize.SpinReserveStates(i)) = 0; %removes state associated with this spinning reserve for this generator.
                        QP.Organize.SpinReserveStates(i) = 0;
                        s_rm = s_rm+1;
                    end

                    if QP.Organize.Equalities(i,1)>0
                        eq = QP.Organize.Equalities(i,:);
                        for j = 1:1:eq(2)
                            reqkeep(eq(1)+(j-1):QP.Organize.t1Balances:(nS-1)*QP.Organize.t1Balances+eq(1)+(j-1)) = 0;%equality constraints
                        end
                        QP.Organize.Equalities(i,1:2) = [0,0];
                        req_rm = req_rm + eq(2);
                    end

                    if isfield(QP.Organize,'Ramping') && QP.Organize.Ramping(i)>0%multi-time step
                        rkeep(QP.Organize.Ramping(i):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.Ramping(i)) = 0;% ramp up constraint
                        rkeep(QP.Organize.Ramping(i)+1:QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.Ramping(i)+1) = 0;% ramp down constraint
                        r_rm = r_rm + 2;
                    end

                    if isfield(QP.Organize,'Inequalities') && QP.Organize.Inequalities(i,1)>0
                        ineq = QP.Organize.Inequalities(i,:);
                        for j = 1:1:ineq(2)
                            rkeep(ineq(1)+(j-1):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+ineq(1)+(j-1)) = 0;%inequality constraints
                        end
                        r_rm = r_rm + ineq(2);
                    end

                    if QP.Organize.SpinRow(i)>0
                        rkeep(QP.Organize.SpinRow(i):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.SpinRow(i)) = 0;% spinning reserve inequality constraint 1
                        rkeep(QP.Organize.SpinRow(i)+1:QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.SpinRow(i)+1) = 0;% spinning reserve inequality constraint 2
                        r_rm = r_rm + 2;
                    end
                else
                    %don't remove, but change state #'s
                    QP.Organize.States(i) = {QP.Organize.States{i}-s_rm}; %empty
                    if QP.Organize.SpinReserveStates(i)>0
                        QP.Organize.SpinReserveStates(i) = QP.Organize.SpinReserveStates(i) - s_rm;
                    end
                    if ~isempty(QP.organize{2,i})
                        QP.organize{2,i} = QP.organize{2,i} -s_rm; %
                    end
                    if i<=nG
                        if QP.Organize.Equalities(i,1)>0
                            QP.Organize.Equalities(i,1) = QP.Organize.Equalities(i,1) - req_rm;
                        end
                        if QP.Organize.Ramping(i)>0%multi-time step
                            QP.Organize.Ramping(i) = QP.Organize.Ramping(i) - r_rm;
                        end

                        if isfield(QP.Organize,'Inequalities') && QP.Organize.Inequalities(i,1)>0
                            QP.Organize.Inequalities(i,1) = QP.Organize.Inequalities(i,1) - r_rm;
                        end

                        if  QP.Organize.SpinRow(i)>0
                            QP.Organize.SpinRow(i) = QP.Organize.SpinRow(i) - r_rm;
                        end
                    end
                end
            end
        end
        QP.Organize.t1States = QP.Organize.t1States - (s_rm - ic_rm);
        QP.Organize.t1ineq = QP.Organize.t1ineq - r_rm;
        QP.Organize.t1Balances = QP.Organize.t1Balances - (req_rm - ic_rm);
        for i = 1:1:n
            if QP.Organize.IC(i)>0
                QP.organize(1,i) = {QP.Organize.IC(i)};
            else 
                QP.organize(1,i) = {[]};%remove disabled generator
            end
            for t= 2:1:nS
                QP.organize(t+1,i) = {QP.organize{t,i}+QP.Organize.t1States};
            end
        end
    end
    for i = nG+1:1:nG+nL
        if  QP.Organize.Transmission(i-nG)>0
            QP.Organize.Transmission(i-nG) = QP.Organize.Transmission(i-nG) - r_rm;
        end
    end
    if isfield(QP.Organize,'HeatVented')
        QP.Organize.HeatVented = max(0,QP.Organize.HeatVented - s_rm);%heat vented states are after all generator & spin reserve states
        %what to do when all generators producing heat at this node are
        %off-line? Remove state?
    end  
end
%% now go through 1 step at a time and remove unnecessary states
% if ~isempty(Locked)     
%     %% I want to re-do this to actually remove states instead of setting lb. Will need to adjust the constraint eqns.
%     for t = 1:1:nS
%         for i = 1:1:nG
%             if Locked(t,i)>0 && Locked(t+1,i) == 0% shutting down
%                 %remove states, and edit ramp down constraint at
%                 %time (t-T) to ramp to a constant, e.g. zero, that
%                 %would signify it shutting down to zero by time t.
%                 %T is how many steps prior it would need to start
%                 %shutting down. Remove all states after t-T.
%             elseif Locked(t,i)==0 && Locked(t+1,i) == 0% shutting down
%                 %remove states, and edit ramp up constraint at
%                 %time (t+T) to ramp to a constant, e.g. zero, that
%                 %would signify it shutting down to zero by time t.
%                 %T is how many steps after now it would reach the 
%                 %lower operating bound. Remove all states until t+T.
%             elseif Locked(t,i)==0
%                 %already shut down, remove states
%             end
%         end
%     end
% end

%% edit QP matrices
rkeep = nonzeros(rkeep);
reqkeep = nonzeros(reqkeep);
xkeep = nonzeros(xkeep);
QP.H = QP.H(xkeep,xkeep);
QP.f = QP.f(xkeep);
QP.Aeq = QP.Aeq(reqkeep,xkeep);
QP.beq = QP.beq(reqkeep);
if r>0 && nnz(rkeep)>0
    QP.A = QP.A(rkeep,xkeep);
    QP.b = QP.b(rkeep);
else QP.A = [];
    QP.b = [];
end
QP.lb = QP.lb(xkeep);
QP.ub = QP.ub(xkeep);

%% old version just lowers upper bound when locking
if ~isempty(Locked)
    index = (0:QP.Organize.t1States:(nS-1)*QP.Organize.t1States)';
    for j = 1:1:nG
        if any(~Locked(:,j)) && any(Locked(:,j)) %is generator off at any point, but not always off and eliminated previously
            Pmin = QP.lb(QP.organize{2,j});
            for t = 1:1:nS
                if ~Locked(t+1,j)
                    QP.lb(QP.organize{t+1,j}) = 0;%was already zero
                    QP.ub(QP.organize{t+1,j}) = 1e-5;%need to make into an equality constraint (first change Aeq(:,i) = 0 to remove from energy balance, if more than 2 concurrent locked off steps, remove ramping constraints)
                    %do this by setting A(:,i)=0 then remove all rows that are all zeros, then rows with 1 value (half a ramping constraint get  a) eliminated and b) converted to an upper bound
                    if isfield(QP.Organize,'SpinReserve')
                        QP.A((t-1)*QP.Organize.t1ineq+QP.Organize.SpinReserve(j),(t-1)*QP.Organize.t1States+QP.Organize.SpinReserveStates(j)) = 0;%%remove from calculation of spinning reserve
                    end
                end
            end
            %% modify lower bounds if ramp rates are to slow (or steps too small)
            rows = QP.Organize.Ramping(j):QP.Organize.t1ineq:(nS-1)*QP.Organize.t1ineq+QP.Organize.Ramping(j);
            states = QP.Organize.States{j};
            nt = length(states); %number of states per timestep
            s = index*ones(1,nt) + ones(nS,1)*states;
            starts = nonzeros((1:nS)'.*((Locked(2:end,j)-Locked(1:nS,j))>0));
            if ~isempty(starts)
                % example, the 4th index of locked is one and it takes 3 steps to turn on (0, .25LB, .7LB, >LB), then starts(1) = 3, n = 3, e = 3. We only need to change LB at 4 & 5 (33% and 67%), 3 was already set to zero
                for k = 1:1:length(starts)
                    e = starts(k);
                    n = 0;
                    Pmax = QP.b(rows(e));
                    for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                        while Pmax<sum(Pmin(1:f)) && (e+n)<=nS
                            QP.lb(s(e+n,f)) = Pmax - sum(QP.lb(s(e+n,1:(f-1)))); %lb to ensure that it is off when locked =0
                            QP.lb(s(e+n,(f+1):nt)) = 0; %other states of this generator during start-up
                            n= n+1; %# of steps to turn on
                            if e+n<=nS
                                Pmax = Pmax+QP.b(rows(e+n));
                            end
                        end
                    end
                end
            end
            stops = nonzeros((1:nS)'.*(Locked(1:nS,j)-(Locked(2:end,j))>0));
            if ~isempty(stops)
                % example, the 6th index of locked is zero and it takes 3 steps to shut down (LB, .67LB, .33LB, 0), then stops(1) = 5, p = 3, e = 3. We only need to change LB at 4 & 5 (67% and 33%), 6 was already set to zero
                for k = 1:1:length(stops)
                    e = stops(k);
                    if e> 1
                        n = 1;
                        Pmax = QP.b(rows(e-1)+1);
                        for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                            while Pmax<sum(Pmin(1:f)) && e-n>0
                                QP.lb(s(e-n,f)) = Pmax - sum(QP.lb(s(e-n,1:(f-1)))); %lb to ensure that it is off when locked =0
                                QP.lb(s(e-n,(f+1):nt)) = 0; %other states of this generator 
                                n= n+1; %# of steps to turn off
                                if e-n>0
                                    Pmax = Pmax+QP.b(rows(e-n)+1);%amount it can ramp down at the previous step
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end 
end%Ends function disableGenerators