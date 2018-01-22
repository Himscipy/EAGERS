function PlotHistogram(h,Data,Xlab)
N = length(Data);
dSort = sort(Data);
Xi = max(1,floor(0.01*N));
Xf = ceil(0.99*N);
range = dSort(Xf)-dSort(Xi);
OoM = log10(range);
if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^(OoM-1);
elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
    Xspace = 10^floor(OoM);
elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
    Xspace = .5*10^floor(OoM);
else  %count in increments of 2, 20, 200 or 2000 etc
    Xspace = .2*10^floor(OoM);
end
hist = [];
label = {};
Xi = ceil(dSort(Xi)/Xspace)*Xspace;
hist(end+1) = nnz(dSort<=Xi);
label(end+1) = {strcat('<',num2str(Xi))};
while Xi<dSort(Xf)
    Xi = Xi + Xspace;
    hist(end+1) = (nnz(dSort<=Xi) - sum(hist(1:end)));
    label(end+1) = {strcat(num2str(Xi-Xspace),'--',num2str(Xi))};
end
hist = hist/N*100;
cla(h)
bar(h,hist)
ylabel(h,'Percent of Time Within Range')
% set(h,'XTickLabel', label)
set(h,'XTickLabel','')
xlim(h,[0.5,length(hist)+.5]);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
pos = get(h,'position');
t = text(0,0,Xlab);
set(t,'Parent',h,'Units','characters','HorizontalAlignment','center','Position',[pos(3)/2,-6,0]);
% Place the text labels
for i = 1:length(hist)
    t = text(0,0,label{i});
    xpos = i*pos(3)/length(hist)- 0.5*pos(3)/length(hist);
    set(t,'Parent',h,'Units','characters','HorizontalAlignment','right','VerticalAlignment','top','Rotation',45,'Position',[xpos,0,0]);
end
end %Ends function PlotHistogram