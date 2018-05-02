function qp = disable_generators(qp,locked,enabled)
% QP are the updated QP matrices
% Organize is the record of which indices are associated with each generator
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Enabled is a variable to include each generator in the optimization or not.
[m,n] = size(qp.organize);
n_s = m-1;
n_g = length(qp.Organize.Dispatchable);
n_b = length(qp.Organize.Building.r);
n_l = length(qp.Organize.Transmission);
if ~isempty(locked)
    enabled = true(n_g,1);
    for j = 1:1:n_g
        if all(locked(:,j)==0)
            enabled(j) = false;
        end
    end
end

%add in constant electric loads from chillers/pumps in FitB
if isfield(qp,'constDemand') 
    Outs = fieldnames(qp.constDemand);
    for s = 1:1:length(Outs)
        if any(qp.constDemand.(Outs{s}).req>0)
            if ~isempty(locked)
                for j = 1:1:n_g
                    qp.constDemand.(Outs{s}).load(:,j) = qp.constDemand.(Outs{s}).load(:,j).*locked(2:end,j);
                end
            elseif ~isempty(enabled)
                qp.constDemand.(Outs{s}).load(:,~enabled) = 0;
            end
            qp.beq(qp.constDemand.(Outs{s}).req) = qp.beq(qp.constDemand.(Outs{s}).req) + sum(qp.constDemand.(Outs{s}).load,2);
        end
    end
end
    
%remove disabled generators from optimization matrices
rmv = false;
for i = 1:1:n_g
    if ~enabled(i)  && ~isempty(qp.Organize.States{i})
        rmv = true;  %there is a disabled generator with states to remove
        break
    end
end

[r,~] = size(qp.A);
[req,x_l] = size(qp.Aeq);
r_keep = linspace(1,r,r)';
req_keep = linspace(1,req,req)';
x_keep = linspace(1,x_l,x_l)';
if rmv    
    %% remove entire generators
    if n_s == 0 %single time step
        s_rm = 0;
        r_rm = 0;
        req_rm = 0;
        for i = 1:1:n
            if ~isempty(qp.Organize.States{i})
                if i<=n_g && ~enabled(i)%remove disabled generator
                    qp.organize{1,i} = [];
                    states = qp.Organize.States{i};
                    x_keep(states) = 0; %removes states associated with this generator.
                    qp.Organize.States(i) = {[]}; %empty
                    s_rm = s_rm+length(states);
                    if qp.Organize.SpinReserveStates(i)>0
                        x_keep(qp.Organize.SpinReserveStates(i)) = 0; %removes state associated with this spinning reserve for this generator.
                        qp.Organize.SpinReserveStates(i) = 0;
                        s_rm = s_rm+1;
                    end
                    
                    if qp.Organize.Equalities(i,1)>0
                        eq = qp.Organize.Equalities(i,:);
                        req_keep(eq(1):eq(1)+eq(2)-1) = 0;%equality constraints
                        qp.Organize.Equalities(i,1:2) = [0,0];
                        req_rm = req_rm + eq(2);
                    end
                    
                    if qp.Organize.SpinRow(i)>0
                        req_keep(qp.Organize.SpinRow(i):qp.Organize.SpinRow(i)+1) = 0;% spinning reserve equality constraint 1 & 2
                        req_rm = req_rm + 2;
                    end
                    if qp.Organize.Inequalities(i)>0
                        r_keep(qp.Organize.Inequalities(i)) = 0;%1 inequality constraint for storage charging state
                        r_rm = r_rm + 1;
                    end
                        
                else
                    %don't remove, but change state #'s
                    qp.Organize.States(i) = {qp.Organize.States{i}-s_rm}; %empty
                    if i<=length(qp.Organize.SpinReserveStates)
                        if qp.Organize.SpinReserveStates(i)>0
                            qp.Organize.SpinReserveStates(i) = qp.Organize.SpinReserveStates(i) - s_rm;
                        end
                    end
                    if ~isempty(qp.organize{1,i})
                        qp.organize{1,i} = qp.organize{1,i} -s_rm; %
                    end
                    if i<=n_g
                        if qp.Organize.Equalities(i,1)>0
                            qp.Organize.Equalities(i,1) = qp.Organize.Equalities(i,1) - req_rm;
                        end
                        if qp.Organize.SpinRow(i)>0
                            qp.Organize.SpinRow(i) = qp.Organize.SpinRow(i) - req_rm;
                        end
                        if qp.Organize.Inequalities(i)>0
                            qp.Organize.Inequalities(i) = qp.Organize.Inequalities(i) - r_rm;
                        end
                    end
                end
            end
        end
        qp.Organize.t1States = qp.Organize.t1States - s_rm;
        qp.Organize.t1ineq = qp.Organize.t1ineq - r_rm;
        qp.Organize.t1Balances = qp.Organize.t1Balances - req_rm;
    else %multi-time step
        ic_rm = 0;%initial conditions removed
        if isfield(qp.Organize,'IC')
            for i = 1:1:n
                if qp.Organize.IC(i)>0
                    if i<=n_g && ~enabled(i)% || (~isempty(Locked) && Locked(1,i)==0))
                        x_keep(qp.Organize.IC(i)) = 0;
                        req_keep(qp.Organize.IC(i)) = 0; %initial condition constraints
                        qp.Organize.IC(i) = 0;
                        ic_rm = ic_rm+1;
                    else
                        qp.Organize.IC(i) = qp.Organize.IC(i)-ic_rm;
                    end
                end
            end
        end
        s_rm = ic_rm;
        r_rm = 0;
        req_rm = ic_rm;
        for i = 1:1:n
            if ~isempty(qp.Organize.States{i})
                if i<=n_g && ~enabled(i)%remove disabled generator
                    qp.organize{2,i} = [];

                    states = qp.Organize.States{i};
                    for j = 1:1:length(states)
                        x_keep(states(j):qp.Organize.t1States:(n_s-1)*qp.Organize.t1States+states(j)) = 0; %removes states associated with this generator.
                    end
                    qp.Organize.States(i) = {[]}; %empty
                    s_rm = s_rm+length(states);

                    if qp.Organize.SpinReserveStates(i)>0
                        x_keep(qp.Organize.SpinReserveStates(i):qp.Organize.t1States:(n_s-1)*qp.Organize.t1States+qp.Organize.SpinReserveStates(i)) = 0; %removes state associated with this spinning reserve for this generator.
                        qp.Organize.SpinReserveStates(i) = 0;
                        s_rm = s_rm+1;
                    end

                    if qp.Organize.Equalities(i,1)>0
                        eq = qp.Organize.Equalities(i,:);
                        for j = 1:1:eq(2)
                            req_keep(eq(1)+(j-1):qp.Organize.t1Balances:(n_s-1)*qp.Organize.t1Balances+eq(1)+(j-1)) = 0;%equality constraints
                        end
                        qp.Organize.Equalities(i,1:2) = [0,0];
                        req_rm = req_rm + eq(2);
                    end

                    if isfield(qp.Organize,'Ramping') && qp.Organize.Ramping(i)>0%multi-time step
                        r_keep(qp.Organize.Ramping(i):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.Ramping(i)) = 0;% ramp up constraint
                        r_keep(qp.Organize.Ramping(i)+1:qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.Ramping(i)+1) = 0;% ramp down constraint
                        r_rm = r_rm + 2;
                    end

                    if isfield(qp.Organize,'Inequalities') && qp.Organize.Inequalities(i,1)>0
                        ineq = qp.Organize.Inequalities(i,:);
                        for j = 1:1:ineq(2)
                            r_keep(ineq(1)+(j-1):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+ineq(1)+(j-1)) = 0;%inequality constraints
                        end
                        r_rm = r_rm + ineq(2);
                    end

                    if qp.Organize.SpinRow(i)>0
                        r_keep(qp.Organize.SpinRow(i):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.SpinRow(i)) = 0;% spinning reserve inequality constraint 1
                        r_keep(qp.Organize.SpinRow(i)+1:qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.SpinRow(i)+1) = 0;% spinning reserve inequality constraint 2
                        r_rm = r_rm + 2;
                    end
                else
                    %don't remove, but change state #'s
                    qp.Organize.States(i) = {qp.Organize.States{i}-s_rm}; %empty
                    if qp.Organize.SpinReserveStates(i)>0
                        qp.Organize.SpinReserveStates(i) = qp.Organize.SpinReserveStates(i) - s_rm;
                    end
                    if ~isempty(qp.organize{2,i})
                        qp.organize{2,i} = qp.organize{2,i} -s_rm; %
                    end
                    if i<=n_g
                        if qp.Organize.Equalities(i,1)>0
                            qp.Organize.Equalities(i,1) = qp.Organize.Equalities(i,1) - req_rm;
                        end
                        if qp.Organize.Ramping(i)>0%multi-time step
                            qp.Organize.Ramping(i) = qp.Organize.Ramping(i) - r_rm;
                        end

                        if isfield(qp.Organize,'Inequalities') && qp.Organize.Inequalities(i,1)>0
                            qp.Organize.Inequalities(i,1) = qp.Organize.Inequalities(i,1) - r_rm;
                        end

                        if  qp.Organize.SpinRow(i)>0
                            qp.Organize.SpinRow(i) = qp.Organize.SpinRow(i) - r_rm;
                        end
                    end
                end
            end
        end
        qp.Organize.t1States = qp.Organize.t1States - (s_rm - ic_rm);
        qp.Organize.t1ineq = qp.Organize.t1ineq - r_rm;
        qp.Organize.t1Balances = qp.Organize.t1Balances - (req_rm - ic_rm);
        for i = 1:1:n
            if qp.Organize.IC(i)>0
                qp.organize(1,i) = {qp.Organize.IC(i)};
            else 
                qp.organize(1,i) = {[]};%remove disabled generator
            end
            for t= 2:1:n_s
                qp.organize(t+1,i) = {qp.organize{t,i}+qp.Organize.t1States};
            end
        end
    end
    for i = n_g+1:1:n_g+n_l
        if  qp.Organize.Transmission(i-n_g)>0
            qp.Organize.Transmission(i-n_g) = qp.Organize.Transmission(i-n_g) - r_rm;
        end
    end
    if isfield(qp.Organize,'HeatVented')
        qp.Organize.HeatVented = max(0,qp.Organize.HeatVented - s_rm);%heat vented states are after all generator & spin reserve states
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
r_keep = nonzeros(r_keep);
req_keep = nonzeros(req_keep);
x_keep = nonzeros(x_keep);
qp.H = qp.H(x_keep,x_keep);
qp.f = qp.f(x_keep);
qp.Aeq = qp.Aeq(req_keep,x_keep);
qp.beq = qp.beq(req_keep);
if r>0 && nnz(r_keep)>0
    qp.A = qp.A(r_keep,x_keep);
    qp.b = qp.b(r_keep);
else qp.A = [];
    qp.b = [];
end
qp.lb = qp.lb(x_keep);
qp.ub = qp.ub(x_keep);

%% old version just lowers upper bound when locking
if ~isempty(locked)
    index = (0:qp.Organize.t1States:(n_s-1)*qp.Organize.t1States)';
    for j = 1:1:n_g
        if any(~locked(:,j)) && any(locked(:,j)) %is generator off at any point, but not always off and eliminated previously
            Pmin = qp.lb(qp.organize{2,j});
            for t = 1:1:n_s
                if ~locked(t+1,j)
                    qp.lb(qp.organize{t+1,j}) = 0;%was already zero
                    qp.ub(qp.organize{t+1,j}) = 1e-5;%need to make into an equality constraint (first change Aeq(:,i) = 0 to remove from energy balance, if more than 2 concurrent locked off steps, remove ramping constraints)
                    %do this by setting A(:,i)=0 then remove all rows that are all zeros, then rows with 1 value (half a ramping constraint get  a) eliminated and b) converted to an upper bound
                    if isfield(qp.Organize,'SpinReserve')
                        qp.A((t-1)*qp.Organize.t1ineq+qp.Organize.SpinReserve(j),(t-1)*qp.Organize.t1States+qp.Organize.SpinReserveStates(j)) = 0;%%remove from calculation of spinning reserve
                    end
                end
            end
            %% modify lower bounds if ramp rates are to slow (or steps too small)
            rows = qp.Organize.Ramping(j):qp.Organize.t1ineq:(n_s-1)*qp.Organize.t1ineq+qp.Organize.Ramping(j);
            states = qp.Organize.States{j};
            nt = length(states); %number of states per timestep
            s = index*ones(1,nt) + ones(n_s,1)*states;
            starts = nonzeros((1:n_s)'.*((locked(2:end,j)-locked(1:n_s,j))>0));
            if ~isempty(starts)
                % example, the 4th index of locked is one and it takes 3 steps to turn on (0, .25LB, .7LB, >LB), then starts(1) = 3, n = 3, e = 3. We only need to change LB at 4 & 5 (33% and 67%), 3 was already set to zero
                for k = 1:1:length(starts)
                    e = starts(k);
                    n = 0;
                    Pmax = qp.b(rows(e));
                    for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                        while Pmax<sum(Pmin(1:f)) && (e+n)<=n_s
                            qp.lb(s(e+n,f)) = Pmax - sum(qp.lb(s(e+n,1:(f-1)))); %lb to ensure that it is off when locked =0
                            qp.lb(s(e+n,(f+1):nt)) = 0; %other states of this generator during start-up
                            n= n+1; %# of steps to turn on
                            if e+n<=n_s
                                Pmax = Pmax+qp.b(rows(e+n));
                            end
                        end
                    end
                end
            end
            stops = nonzeros((1:n_s)'.*(locked(1:n_s,j)-(locked(2:end,j))>0));
            if ~isempty(stops)
                % example, the 6th index of locked is zero and it takes 3 steps to shut down (LB, .67LB, .33LB, 0), then stops(1) = 5, p = 3, e = 3. We only need to change LB at 4 & 5 (67% and 33%), 6 was already set to zero
                for k = 1:1:length(stops)
                    e = stops(k);
                    if e> 1
                        n = 1;
                        Pmax = qp.b(rows(e-1)+1);
                        for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                            while Pmax<sum(Pmin(1:f)) && e-n>0
                                qp.lb(s(e-n,f)) = Pmax - sum(qp.lb(s(e-n,1:(f-1)))); %lb to ensure that it is off when locked =0
                                qp.lb(s(e-n,(f+1):nt)) = 0; %other states of this generator 
                                n= n+1; %# of steps to turn off
                                if e-n>0
                                    Pmax = Pmax+qp.b(rows(e-n)+1);%amount it can ramp down at the previous step
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end 
end%Ends function disable_generators