function Profile = load_sched2(Data,Date)
%Interpolates two fields Data.Timestamp and Data.Load
%If the date is in the timestamp range, it uses the value
%If not it tries to use the same day of the 1st year in the timestamp
%If that day doesn't exist for a different year, it maches the day of the month,
D = datevec(Date(1));
dt = Date(2)-Date(1);
days = max(1,ceil(Date(end) - datenum([D(1),D(2),D(3)])));
nS = length(Date); % number timesteps
Profile = zeros(nS,1);
[newDay,D_of_Y, D_of_M] = day_of_year(Date);
r = (Data.Timestamp(2)-Data.Timestamp(1))/(Date(2)-Date(1));% # of points that must be created between datum 
for d = 1:1:days
    p = newDay(d);
    if d<days
        p2 = newDay(d+1)-1;
    else 
        p2 = length(Date);
    end
    if Date(p)>=floor(Data.Timestamp(1)) && (Date(p2)-dt)<floor(Data.Timestamp(end))
        %Date is within the data period, use directly
        if Data.Timestamp(1)>(Date(p)-dt)
            Xi = 1;
        else
            Xi = nnz(Data.Timestamp<=Date(p)-dt);
        end
        Xf = nnz(Data.Timestamp<=(Date(p)+1));
    else%find matching day in different year
        [data_newDay,data_D_of_Y,data_D_of_M] = day_of_year(Data.Timestamp);
        data_days = length(data_newDay);
        if any(data_D_of_Y == D_of_Y(d))
            Index = nonzeros((1:data_days)'.*(data_D_of_Y == D_of_Y(d)));
            Xi = data_newDay(Index(1));
            if Index(1) == data_days
                Xf = length(Data.Timestamp);
            else
                Xf = data_newDay(Index(1)+1)-1;
            end
        elseif any(data_D_of_M == D_of_M(d))
            Index = nonzeros((1:data_days)'.*(data_D_of_M == D_of_M(d)));
            Xi = data_newDay(Index(1));
            if Index(1) == data_days
                Xf = length(Data.Timestamp);
            else
                Xf = data_newDay(Index(1)+1)-1;
            end
        else
            rn = ceil(rand(1)*data_days);
            Xi = data_newDay(rn);
            if rn == data_days
                Xf = length(Data.Timestamp);
            else
                Xf = data_newDay(rn+1)-1;
            end
        end
    end
    
    if abs(r-1)<1e-8
        Profile(p:p2,1) = Data.Load(Xi:Xf);
    elseif r<1 %extra datum, average points in between
        nS2 = length(Data.Timestamp);
        x1 = Xi;
        for k = p:1:p2
            if k == p2
                x2 = Xf;
            else
                x2 = x1;
                while x2<nS2 && (Data.Timestamp(x2+1)-floor(Data.Timestamp(x2+1)))<(Date(k)-floor(Date(k)))
                    x2 = x2+1;
                end
            end
            Profile(k,1) = mean(Data.Load(x1:x2));
            x1 = x2+1;
        end
    elseif r>1 %interpolate between timesteps
        t = Data.Timestamp(Xi:Xf) - floor(Data.Timestamp(Xi:Xf));
        t2 = Date(p:p2)-floor(Date(p));
        Profile(p:p2,1) = interp1(t,Data.Load(Xi:Xf),t2);
    end
end
end%Ends function load_sched2

function [newDay,D_of_Y, D_of_M] = day_of_year(Date)
D = datevec(Date(1));
sd = floor(Date(1)); % start date
days = max(1,ceil(Date(end) - datenum([D(1),D(2),D(3)])));
newDay = zeros(days,1);
D_of_Y = zeros(days,1);
D_of_M = zeros(days,1);
p = 1;
nS = length(Date);
for d = 1:1:days
    p2 = p;
    while p2<nS && Date(p2+1)<=sd+d
        p2 = p2+1;
    end
    newDay(d) = p;
    D = datevec(Date(p));
    D_of_Y(d) = floor(Date(p) - datenum([D(1),1,1]))+1;
    D_of_M(d) = floor(Date(p) - datenum([D(1),D(2),1]))+1;
    p = p2+1;
end

end%Ends function day_of_year