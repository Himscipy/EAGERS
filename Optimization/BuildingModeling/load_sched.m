function Out = load_sched(Build,Date,param)
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
        Sched = convert_sched(Build.Schedule.(param).(day){s},Ramp);
    else
        Sched = convert_sched(Build.Schedule.(param).(day),Ramp);
    end
    Out(p:p2,1) = interp1(Sched(:,1),Sched(:,2),hour((p:p2)));
    p = p2+1;
    wd = wd+1;
end
end%Ends function load_sched