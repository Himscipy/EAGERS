function PlotCosts(handles)
global testSystems
if isfield(handles,'LegendDeleteProxy')%2013 matlab
    delete(handles.LegendColorbarLayout)
    delete(handles.LegendDeleteProxy)
elseif isfield(handles,'legend')%2015 matlab
    delete(handles.legend)
end
h1 = handles.axesMain;
cla(h1);
h2 = handles.axesCumulative;
cla(h2);
nShow = length(testSystems); %number of systems to plot
Timestamp = testSystems(1).Design.Timestamp;
labelMonths = {};
m = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';};
D = datevec(Timestamp(1));
n = 0;
while datenum([D(1) D(2)+n 1])<Timestamp(end)
    labelMonths(end+1) = m(rem(D(2)+n-1,12)+1);
    n = n+1;
end
Xpos = zeros(n,nShow);
if nShow == 1
    Xpos(:,1) = (.5:1:n)';
else
    for i = 1:1:nShow
        Xpos(:,i) = ((1/nShow^2 + (i-1)/nShow):1:n)';
    end
end
NPC = zeros(nShow,1);
Name = cell(nShow,1);
colormap(h1,'summer')
for i = 1:1:length(testSystems)
    Data = testSystems(i).Costs.Design;
    NPC(i) = testSystems(i).Costs.NPC;
    Name(i) = {testSystems(i).Name};
    %savebuilding color and figure out spacing/width
    bar(h1,Xpos(:,i),Data,'stacked','BarWidth',0.8/nShow)
%     set(h1,'ColorOrderIndex',1);
end
legend(h1,{'Financing Charges';'O & M Charges';'Re-start Charges';'Demand Charges';'Electric Use Charges';'Fuel Charges';})
xlim(h1,[0,n])
ylabel(h1,'Cost ($)','Color','k','FontSize',14)
set(h1,'XTick',mean(Xpos,2),'XTickLabel', labelMonths)
colormap(h2,'summer')
bar(h2,Xpos(1,:),NPC)
xlim(h2,[0,1])
ylabel(h2,'20 Year Net Present Cost ($)','Color','k','FontSize',14)
set(h2,'XTick',Xpos(1,:),'XTickLabel',Name) 
x_lab = strcat('Order is : ',testSystems(1).Name);
for i = 2:1:nShow
    x_lab = strcat(x_lab,' , ',testSystems(i).Name);
end
xlabel(h1,x_lab)