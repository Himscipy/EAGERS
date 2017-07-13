function ExportResult
global Model_dir Plant
if isfield(Plant,'RunData') && ~isempty(Plant.RunData)
    [f,p]=uiputfile(fullfile(Model_dir,'results','LoggedData','Test1.mat'),'Save logged test data...');
    save([p,f],'Plant.RunData')
end
K = menu('Save Dialog','Export Dispatch to Excel','Do not Export');
nG = length(Plant.Generator);
if K ==1
    [f,p]=uiputfile(fullfile(Model_dir,'results',strcat(Plant.Name,'.xls')),'Save Plant As...');
    if f==0; return; end
    filename = fullfile(p,f);
    MSG = msgbox('Exporting...');

    ColumnLabelSheet1 = {'Name';'Type';'Size';};
    xlswrite(filename,ColumnLabelSheet1,'System','A1:A3');
    Row1Sheet1 = {Plant.Name;};
    Row2Sheet1 = {'Fuel  -->'};
    Row3Sheet1 = {'Size (kW) -->'};
    Row4Sheet1 = {'Day'};

    D = datevec(Plant.Dispatch.Timestamp');
    Data = D(:,3);
    Row4Sheet1(end+1) = {'Hour'};
    Data(:,end+1) = D(:,4) + round(4*D(:,5)/60)/4;
    Row4Sheet1(end+1) = {'Temperature'};
    Data(:,end+1) = Plant.Dispatch.Temperature;
    if isfield(Plant.Dispatch.Demand,'E')
        Row4Sheet1(end+1) = {'Electric Demand (kW)'};
        Data(:,end+1) = Plant.Dispatch.Demand.E';
    end
    if isfield(Plant.Dispatch.Demand,'H')
        Row4Sheet1(end+1) = {'Heating Demand (kW)'};
        Data(:,end+1) = Plant.Dispatch.Demand.H';
    end
    if isfield(Plant.Dispatch.Demand,'C')
        Row4Sheet1(end+1) = {'Cooling Demand (kW)'};
        Data(:,end+1) = Plant.Dispatch.Demand.C';
    end
    s = size(Data);
    j = s(2);
    for i = 1:1:nG
        j = j+1;
        Row1Sheet1(j) = Plant.Generator(i).Name;
        Row2Sheet1(j) = Plant.Generator(i).Source;
        Row3Sheet1(j) = Plant.Generator(i).Size;
        Data(:,j) = Plant.Dispatch.GeneratorState(:,i);
        if strfind(Plant.Generator(i).Type,'Utility')
            Row4Sheet1(j) = {'Energy Purchase (kW)'};
        elseif strfind(Plant.Generator(i).Type,'Storage')
            Row4Sheet1(j) = {'Usable Energy Stored (kWh)'};
        else %generators, chiller, heaters
            Row4Sheet1(j) = {'Ouput (kW)'};
            j = j+1;
            Row4Sheet1(j) = {'Input (kW)'};
            Data(:,j) = Plant.Dispatch.GeneratorInput(:,i);
        end
    end
    xlswrite(filename,Row1Sheet1,'System','B1');
    xlswrite(filename,Row2Sheet1,'System','B2');
    xlswrite(filename,Row3Sheet1,'System','B3');
    xlswrite(filename,Row4Sheet1,'System','B4');
    [m,n] = size(Data);
    Data2 = {};
    for i = 1:1:n
        Data2(:,i) = cellstr(num2str(Data(:,i),3));
    end
    xlswrite(filename,Data2,'System','B5');
    close(MSG)
end