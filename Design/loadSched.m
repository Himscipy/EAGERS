function Out = loadSched(Build,Date,param)
D = datevec(Date(1));
days = max(1,ceil(Date(end) - datenum([D(1),D(2),D(3)])));
sd = floor(Date(1)); % start date
nS = length(Date); % number timesteps
Out = zeros(nS,1);
[~,~,~,H,M,S] = datevec(Date);
H(H==0&M==0&S==0) = 24; % change hour 0 to hour 24 for interpolating.
hour = H+M/60+S/3600;
daysAfter1_1_17 = floor(Date(1)) - datenum([2017,1,1]);
wd = 1 + mod(daysAfter1_1_17,7); % day of the week, sunday is 1, saturday is 7
p = 1;
for d = 1:1:days
    p2 = p;
    while p2<nS && Date(p2+1)<=sd+d
        p2 = p2+1;
    end
    D = datevec(Date(p));
    h_of_y = 24*(Date(p) - datenum([D(1),1,1])); % hour since start of year
    s = nnz(Build.Schedule.(param).Seasons<=h_of_y)+1;% get season
    if isfield(Build.Schedule.(param),'Ramp')
        Ramp = Build.Schedule.(param).Ramp;
    else
        Ramp = 1e-4;
    end
    % get schedule
    if wd == 1% Sunday
        day = 'Sun';
    elseif wd == 7 % Saturday
        day = 'Sat';
        wd = wd-7;
    else % Weekday
        day = 'Weekday';
    end
    if iscell(Build.Schedule.(param).(day))
        Sched = convertSched(Build.Schedule.(param).(day){s},Ramp);
    else
        Sched = convertSched(Build.Schedule.(param).(day),Ramp);
    end
    Out(p:p2,1) = interp1(Sched(:,1),Sched(:,2),hour((p:p2)));
    p = p2+1;
    wd = wd+1;
end
end%Ends function loadSched


function newsched = convertSched(sched,Ramp)
[m,n] = size(sched);
if m==2
    newsched = sched; %this is constant all day, already made 0 hour and 24 hour
else
    newsched = zeros(2*(m-1),n);
    sched(1,2) = sched(2,2);
    newsched(1,:) = sched(1,:);%hour 0
    newsched(end,:) = sched(m,:);%hour 24
    if Ramp<1e-3
        for i = 2:1:m-1 
            newsched(2*i-2,1) = sched(i);
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+Ramp;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    else 
        for i = 2:1:m-1 %add points in the middle so it can be properly interpolated
            t_bef = sched(i) - sched(i-1);
            t_aft = sched(i+1)-sched(i);
            Ramp2 = min([Ramp, t_aft/2,t_bef/2]);
            newsched(2*i-2,1) = sched(i)-0.5*Ramp2;
            newsched(2*i-2,2) = sched(i,2);
            newsched(2*i-1,1) = sched(i)+0.5*Ramp2;
            newsched(2*i-1,2) = sched(i+1,2);
        end
    end
end
end%Ends function convertSched