
DamNames = {'Grand Coulee';'Chief Joseph';'Wells';'Rocky Reach';'Rock Island';...
            'Wanapum';'Priest Rapids';'McNary';'John Day';'The Dalles';'Bonneville';...
            'Lower Granite';'Little Goose';'Lower Monumental';'Ice Harbor'};
        
NodeNames = {'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa';'Wenatchee, Wa';...
             'South Wenatchee, Wa';'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';...
             'Rufus, Or';'The Dalles, Or';'Bonneville, Or';'Almota, Wa'; 'Starbuck, Wa';...
             'Kahlotus, Wa';'Pasco, Wa'};

UpRiver = {'';'Grand Coulee, Wa';'Bridgeport, Wa';'Azwell, Wa';'Wenatchee, Wa';...
            'South Wenatchee, Wa';...
             'Beverly, Wa';{'Mattawa, Wa';'Pasco, Wa'};'Umatilla, Or';'Rufus, Or';...
             'The Dalles, Or';'Bonneville, Or'; '';'Starbuck, Wa';...
             'Kahlotus, Wa';'Pasco, Wa';'Umatilla, Or'};
         
DownRiver = {'Bridgeport, Wa';...
             'Azwell, Wa';'Wenatchee, Wa';'South Wenatchee, Wa';...
             'Beverly, Wa';'Mattawa, Wa';'Umatilla, Or';'Rufus, Or';...
             'The Dalles, Or';'Bonneville, Or'; '';'Starbuck, Wa';...
             'Kahlotus, Wa';'Pasco, Wa';'Umatilla, Or'};

%% Loading Files & ReOrginizing
%These are corrected data 
load('SourceSinkHrly')
load('InFlowCrctd')
load('OutFlowNMD')
load('powerFlowNMD')

for i = 1:1:length(DamNames) 
    if exist(strcat(DamNames{i},'_','2007','_','2016','.mat'),'file')
        name = strcat(DamNames{i},'_','2007','_','2016');
        R = load(name);
        ix = fieldnames(R);
        Dams(i,1) = { R.(ix{1}) };
    else Dams{i,1} = {};
    end
end 


%% Calculating Inflow-Outflow

for i = 1:1:length(DamNames)
    InOut(:,i) = crctdInFlow(:,i)-OutFlow(:,i);
    if i == 1 || i == 12
        InOut(:,i) = SourceSinkHrly(:,i)-OutFlow(:,i);
    end
end

%% Calculating Yearly totals

year = {'2007';'2008';'2009';'2010';'2011';'2012';'2013';'2014';'2015';'2016'};
month = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'};
normDays = [31;28;31;30;31;30;31;31;30;31;30;31];
leapDays = [31;29;31;30;31;30;31;31;30;31;30;31];

for i = 1:1:length(DamNames)
    x = 0;
    for y = 1:1:length(year)
        for m = 1:1:length(month)
            if strcmp(year{y},'2008') || strcmp(year{y},'2012') || strcmp(year{y},'2016') 
                hours = leapDays(m)*24;
            else
                hours = normDays(m)*24;
            end
           monthly{i}(m,y) = sum(InOut(x+1:x+hours,i));
           x = x + hours;
        end 
    end
end
