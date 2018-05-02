function PlotCosts(handles,names,timestamp,costs,npc)
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
n_sys = length(names); %number of systems to plot
labelMonths = {};
m = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';};
D = datevec(timestamp(1));
n = 0;
while datenum([D(1) D(2)+n 1])<timestamp(end)
    labelMonths(end+1) = m(rem(D(2)+n-1,12)+1);
    n = n+1;
end
Xpos = zeros(n,n_sys);
if n_sys == 1
    Xpos(:,1) = (.5:1:n)';
else
    for i = 1:1:n_sys
        Xpos(:,i) = ((1/n_sys^2 + (i-1)/n_sys):1:n)';
    end
end
colormap(h1,'summer')
for i = 1:1:n_sys
    Data = costs(:,:,i);
    %savebuilding color and figure out spacing/width
    bar(h1,Xpos(:,i),Data,'stacked','BarWidth',0.8/n_sys)
    set(h1,'ColorOrderIndex',1);
end
legend(h1,{'Financing Charges';'O & M Charges';'Re-start Charges';'Demand Charges';'Electric Use Charges';'Fuel Charges';})
xlim(h1,[0,n])
ylabel(h1,'Cost ($)','Color','k','FontSize',14)
set(h1,'XTick',mean(Xpos,2),'XTickLabel', labelMonths)
colormap(h2,'summer')
bar(h2,Xpos(1,:),npc)
xlim(h2,[0,1])
ylabel(h2,'20 Year Net Present Cost ($)','Color','k','FontSize',14)
set(h2,'XTick',Xpos(1,:),'XTickLabel',names) 
x_lab = strcat('Order is : ',names{1});
for i = 2:1:n_sys
    x_lab = strcat(x_lab,' , ',names{i});
end
xlabel(h1,x_lab)