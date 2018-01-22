function Data = GetCurrentData(Date)
global RealTimeData RealTime Virtual
if ~isinf(Date(1)) &&  all(Date<=RealTimeData.Timestamp(end))
    %% need better indexing system to call current data
    F = fieldnames(RealTimeData);
    F = F(~strcmp('Timestamp',F));
    F = F(~strcmp('Holidays',F));
    x1 = max(1,nnz(RealTimeData.Timestamp<Date(1))); 
    Data.Timestamp = Date;
    if length(Date) == 1
        r = (Date - RealTimeData.Timestamp(x1))/(RealTimeData.Timestamp(x1+1) - RealTimeData.Timestamp(x1));
        for j = 1:1:length(F)
            if isstruct(RealTimeData.(F{j}))
                S = fieldnames(RealTimeData.(F{j}));
                for i = 1:1:length(S)
                    Data.(F{j}).(S{i}) = (1-r)*RealTimeData.(F{j}).(S{i})(x1,:) + r*RealTimeData.(F{j}).(S{i})(x1+1,:);
                end
            else
                Data.(F{j}) = (1-r)*RealTimeData.(F{j})(x1) + r*RealTimeData.(F{j})(x1+1);
            end
        end        
    else
        x2 = max(1,nnz(RealTimeData.Timestamp<Date(end))); 
        for j = 1:1:length(F)
            if isstruct(RealTimeData.(F{j}))
                S = fieldnames(RealTimeData.(F{j}));
                for i = 1:1:length(S)
                    Data.(F{j}).(S{i}) = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.(F{j}).(S{i})(x1:x2+1,:),Date);
                    Data.(F{j}).(S{i})(abs(Data.(F{j}).(S{i}))<1e-2) = 0;
                end
            else
                Data.(F{j}) = interp1(RealTimeData.Timestamp(x1:x2+1),RealTimeData.(F{j})(x1:x2+1,:),Date);
                Data.(F{j})(abs(Data.(F{j}))<1e-2) = 0;
            end
        end
    end
else % End Operation or Simulation
    if RealTime
        closePorts;
    end
    RealTime=0;%end condition for real simulation
    Virtual = 0;%end condition for virtual simulation
    
    T1 = timerfind('Name', 'dispTimer') ;
    T2 = timerfind('Name', 'optTimer') ;
    T3 = timerfind('Name', 'mpcTimer') ;
    T4 = timerfind('Name', 'fanTimer') ;
    Timers = [T1,T2,T3,T4];
    for i = 1:1:length(Timers)
        stop(Timers(i));
        delete(Timers(i))
    end
    Data = [];
end