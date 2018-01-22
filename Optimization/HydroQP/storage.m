clear all
clear 

Model_dir = strrep(which('EAGERS.m'),'EAGERS.m',''); %Have to run eagers at least once to locate
load(fullfile(Model_dir,'Data','Hydro','SourceSinkHrly.mat'))
load(fullfile(Model_dir,'Data','Hydro','InFlowCrctd.mat'))
load(fullfile(Model_dir,'Data','Hydro','OutFlowNMD.mat'))
load(fullfile(Model_dir,'Data','Hydro','powerFlowNMD.mat'))
% load('Headdata.mat');

folder = 'C:\Users\MME-Admin\Documents\GitHub\EAGERS_Internal\Optimization\HydroQP\Head'; %location to save graphs

DamNames = {'Grand Coulee', 'Chief Joseph', 'Wells', 'Rocky Reach', 'Rock Island', ...
        'Wanapum', 'Priest Rapids', 'McNary', 'John Day', 'The Dalles', ...
        'Bonneville', 'Lower Granite', 'Little Goose', 'Lower Monumental',...
        'Ice Harbor'};
nD = length(DamNames);

for i = 1:1:nD
    if exist(fullfile(Model_dir,'Data','Hydro',strcat(DamNames{i},'_','2007','_','2016','.mat')),'file')
        name = fullfile(Model_dir,'Data','Hydro',strcat(DamNames{i},'_','2007','_','2016'));
        R = load(name);
        ix = fieldnames(R);
        Dams(i,1) = { R.(ix{1}) };
    else Dams{i,1} = {};
    end
end 

for i = 1:1:nD
    Head(:,i) = Dams{i}(:,7);
end

for i = 1:1:nD
    for k = 1:1:length(Head(:,i))
        if isnan(Head(k,i)) 
            Head(k,i) = Head(k-1,i);
        elseif (Head(k,i)) <= 0
            Head(k,i) = Head(k-1,i);
%         elseif (Head(k,i)) >= 100
%             Head(k,i) = Head(k-1,i);
        end
    end
end

Timestamp = Dams{1}(:,1:2);
adjDate = datenum([1899 12 30 1 0 0]); %adjust excel format to matlab
nS = length(Timestamp);
z = 1;
q = 1;
for i = 1:1:nS
    dateFull = datevec(Timestamp(i,1) + adjDate);
    
    %Filling out yearly
    if i == 1
        Year(z,1) = dateFull(1);
        Year(z,2) = i; %start of first year
        z = z+1;
    elseif dateFull(1) > Year(z-1,1)
        Year(z-1,3) = i-1; %end of previous year
        Year(z,1) = dateFull(1);
        Year(z,2) = i; %start of this year
        z = z+1;
    elseif i == nS
        Year(z-1,3) = i; %end of previous year
    end 
    
    %Filling out monthly
    if i == 1     
        Month(q,1) = dateFull(2);
        Month(q,2) = i;
        q = q+1;
    elseif dateFull(2) > Month(q-1,1)
        Month(q-1,3) = i-1;
        Month(q,1) = dateFull(2);
        Month(q,2) = i;
        q = q+1;
    elseif i == nS
        Month(q-1,3) = i;
    end
end 
nYears = length(Year);

loadGen = [378.9256198	21.32231405	13.68595041	15.78512397	5.41322314	32.89256198	9.797520661	55.78512397	104.5454545	13.63636364	22.19008264	18.19008264	21.33471074	17.85123967	10.2892562];
% minsvals = [240; 155; 45; 70; 25; 40; 60; 64; 90; 65; 35; 93; 90; 90; 80];
% maxsvals = [340; 185; 80; 100; 65; 95; 90; 78; 110; 87; 70; 105; 105; 105; 105];
currentDays = (Year(:,3)-Year(:,2))+1;
mostDays = max(currentDays);
for j = 1:1:nD
    HeadYear = zeros(mostDays,nYears);
    for i = 1:1:nYears
        HeadYear(1:currentDays(i,1),i) = Head(Year(i,2):Year(i,3),j);
    end  
    minsvals = min(Head(:,j));
    maxsvals = max(Head(:,j));
    avgvals = mean(Head(:,j));

    storageHalf(1,j) = ((maxsvals-minsvals)/2)+minsvals;
    storageHalf(2,j) = avgvals;
    storageHalf(3,j) = loadGen(1,j);
    storageHalf(4,j) = (mean(Head(:,j))./max(Head(:,j))).*max(Head(:,j));
    
    plot(1:mostDays,HeadYear)
    ylim([minsvals maxsvals]);
    legend('2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','maxMin','Avg','loadGen','Location','southeast');
    xticks([Month(1,2),Month(2,2),Month(3,2),Month(4,2),Month(5,2),Month(6,2),Month(7,2),Month(8,2),Month(9,2),Month(10,2),Month(11,2),Month(12,2)]);
    xticklabels({'Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec'})
    
    hold on
    plot(1:mostDays,ones(mostDays,1)*(storageHalf(1,j)),'<c')
    plot(1:mostDays,ones(mostDays,1)*(storageHalf(2,j)),'oy')
    plot(1:mostDays,ones(mostDays,1)*(storageHalf(3,j)),'*g')
    plot(1:mostDays,ones(mostDays,1)*(storageHalf(4,j)),'--*k')
    hold off
    
    saveas(gcf, fullfile(folder, sprintf('Dam Years Head %d.jpg', j)))
end 

for i = 1:1:nYears
    HeadYear = [];
    HeadYear = zeros(currentDays(i,1),nD);
    for j = 1:1:nD
        HeadYear(1:currentDays(i,1),j) = Head(Year(i,2):Year(i,3),j);
    end  
    plot(1:currentDays(i,1),HeadYear)
    legend(DamNames,'Location','southeast');
    xticks([Month(1,2),Month(2,2),Month(3,2),Month(4,2),Month(5,2),Month(6,2),Month(7,2),Month(8,2),Month(9,2),Month(10,2),Month(11,2),Month(12,2)]);
    xticklabels({'Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec'})
    saveas(gcf, fullfile(folder, sprintf('Yearly Head %d.jpg', Year(i,1))))
end 
