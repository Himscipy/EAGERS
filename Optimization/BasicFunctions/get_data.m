function select_data = get_data(data,date,param)
%load current or past data
n_s = length(date);
need_interp = true;
index = zeros(n_s,1);
select_data.Timestamp = date;
if date(1)<data.Timestamp(1)
    date = date + ceil(data.Timestamp(1)-date(1));%add a whole # of days
end
index(1) = nnz(data.Timestamp<date(1)+1e-6); %+1e-6 to avoid rounding problems
if n_s==1
    need_interp = false;
elseif index(1)>0 && abs(data.Timestamp(index(1))-date(1))<1e-6 && abs(data.Timestamp(index(1)+1)-date(2))<1e-6 %no need to interpolate
    index(2:n_s) = index(1) + (1:n_s-1)';
    need_interp = false;
else
    x1 = index(1);
    x2 = max(1,nnz(data.Timestamp<date(end)));
end

if ~isempty(param)
    if need_interp
        if length(param) == 1
            select_data = interp1(data.Timestamp(x1:x2+1),data.(param{1})(x1:x2+1,:),date);
        elseif length(param) == 2
            select_data = interp1(data.Timestamp(x1:x2+1),data.(param{1}).(param{2})(x1:x2+1,:),date);
        end
    else
        if length(param) == 1
            select_data = data.(param{1})(index,:);
        elseif length(param) == 2
            select_data = data.(param{1}).(param{2})(index,:);
        end
    end
else
    f = fieldnames(data);
    f = f(~strcmp('Timestamp',f));
    for j = 1:1:length(f)
        if isstruct(data.(f{j}))
            s = fieldnames(data.(f{j}));
            for i = 1:1:length(s)
                if isnumeric(data.(f{j}).(s{i})) && ~isempty(data.(f{j}).(s{i}))
                    if need_interp
                        select_data.(f{j}).(s{i}) = interp1(data.Timestamp(x1:x2+1),data.(f{j}).(s{i})(x1:x2+1,:),date);
                    else
                        select_data.(f{j}).(s{i}) = data.(f{j}).(s{i})(index,:);
                    end
                end
            end
        else
            if isnumeric(data.(f{j})) && ~isempty(data.(f{j}))
                if need_interp
                    select_data.(f{j}) = interp1(data.Timestamp(x1:x2+1),data.(f{j})(x1:x2+1,:),date);
                else
                    select_data.(f{j}) = data.(f{j})(index,:);
                end
            end
        end
    end 
end
end%Ends function get_data