%% create graphs for comparison from Matlab
clear all
clc
load('newUsableSize.mat')
% load('30days_distributedDemand_20132014_April.mat')
load('forecast.mat')
% load('30days.mat') %This is 30 days, 4 hr timestamps, all dams/nodes
% load('plantData.mat') %This is full year, 4 hr timestamps, all dams/nodes
% folder = 'C:\Users\MME-Admin\Documents\Figures\185 timesteps Distributed Demand 20132014\April'; %location to save graphs
 folder = 'C:\Users\MME-Admin\Documents\Figures\newUsableSize';
 
nG = length(Plant.Generator); %number of Generators
lineNames = Plant.subNet.Electrical.lineNames;
nL = length(lineNames);
OptData = Plant.Dispatch; %Where optimization data is stored
Historical = Plant.Data.Hydro; %Where historical Data is stored
HDate = Historical.Timestamp;
% Date = [];
% for i = 1:1:length(OptData.Timestamp)-1
ptime = (7*24)/4;
Date = OptData.Timestamp(ptime+1:2*ptime,1);
% for i = 1:1:ptime
%     if nnz(OptData.Timestamp(i,1))
%         Date(end+1,1) = OptData.Timestamp(i,1);
%     end
% end

dt = floor(mean((Date(2:end) - Date(1:end-1))*24));
nS = length(Date);


D = datevec(Date);
HD = datevec(HDate);
for i = 1:1:length(D)
    D(i,4) = D(i,4) + ceil(D(i,5)/60);
end
for i = 1:1:length(HD)
    HD(i,4) = HD(i,4) + ceil(HD(i,5)/60);
end

%Calculate balances for graphing
networkNames = fieldnames(Plant.subNet); 
for net = 1:1:length(networkNames)
    nN = length(Plant.subNet.(networkNames{net}).nodes);
    if strcmp((networkNames{net}),'Electrical') %Electric nodal balance, transfers, and losses    
        %% for Real Historical Electric Data
        xi = nnz(Historical.Timestamp<Date(1))+1; %Start date for actual data
        xf = nnz(Historical.Timestamp<Date(end))+1; %End date for actual data
        r = round((Date(2)-Date(1))/(Historical.Timestamp(2) - Historical.Timestamp(1)));
        H = {'PowerGen';'OutFlow';'InFlow';'PowerFlow';'SpillFlow';'Head';};
        z = 1;
        x1 = xi;
        while x1<=xf
            x2 = x1-1+r;
            for j = 1:1:length(H)
                for i = 1:1:nG
                    Hist.(H{j})(z,i) = mean(Historical.(H{j})(x1:x2,i));
                end
            end
            Hist.Timestamp(z,1) = Historical.Timestamp(x1,1);
            x1 = x2+1;
            z = z+1;
        end
%         for i = 1:1:length(Historical.Timestamp)
%             if (D(1,1) == HD(i,1)) && (D(1,2) == HD(i,2)) && (D(1,3) == HD(i,3)) && (D(1,4) == HD(i,4)) && isempty(xi)
%                 xi = i;
%             elseif (D(end,1) == HD(i,1)) && (D(end,2) == HD(i,2)) && (D(end,3) == HD(i,3)) && (D(end,4) == HD(i,4))
%                 xf = i;
%             end
%         end

%         HistAGen = zeros(nS,nN);
%         Hist.PowerGen = zeros(nS,nG);
% %         HistPFGen = zeros(nS,nG);
%         Hist.OutFlow = zeros(nS,nG);
%         Hist.InFlow = zeros(nS,nG);
%         Hist.PowerFlow = zeros(nS,nG);
%         Hist.SpillFlow = zeros(nS,nG);
%         Hist.Head = zeros(nS,nG);
%         for i = 1:1:nG
%             n = Plant.Generator(i).QPform.(networkNames{net}).subnetNode; %Electric node number
%             
%             for t = xi:4:xf
%                 HistAGen(z,n) = HistAGen(z,n) + mean(Historical.PowerGen(t:t+3,i));%PowerGenData at each node at each timestep
%                 Hist.PowerGen(z,i) = Hist.PowerGen(z,i) + mean(Historical.PowerGen(t:t+3,i));%PowerGenData at each Dam at each timestep
% %                 HistPFGen(z,i) = HistPFGen(z,i) + mean(Historical.PowerFlow(t:t+3,i))/Plant.Generator(i).QPform.Stor.Power2Flow;%PowerGenData at each Dam at each timestep for powerflow conversion
% %                 HistDemand(z,n) = mean(Plant.Data.Demand.E(t:t+3,i)); %mean of stored demand is our demand
%                 Hist.OutFlow(z,i) = Hist.OutFlow(z,i) + mean(Historical.OutFlow(t:t+3,i));%Historical OutFlow data
%                 Hist.InFlow(z,i) = Hist.InFlow(z,i) + mean(Historical.InFlow(t:t+3,i));%Historical InFlow data
%                 Hist.PowerFlow(z,i) = Hist.PowerFlow(z,i) + mean(Historical.PowerFlow(t:t+3,i));%Historical PowerFlow data
%                 Hist.SpillFlow(z,i) = Hist.SpillFlow(z,i) + mean(Historical.SpillFlow(t:t+3,i));%Historical SpillFlow data
%                 Hist.Head(z,i) = Hist.Head(z,i) + mean(Historical.Head(t:t+3,i));%Historical Head data
%                 z = z+1;
%             end
%         end
        
        %% for optimization data
        nodalEnergyBalance = zeros(nS,nN);    
        Generated = zeros(nS,nN);
        Demand = zeros(nS,nN);
        for i = 1:1:nG % Summing the amount of Energy produced at each node 
            n = Plant.Generator(i).QPform.(networkNames{net}).subnetNode; %Electric node number
            Generated(:,n) = Generated(:,n) + OptData.GeneratorState(1:nS,i); %Sum all power Generated at that node
            if ~isempty(Plant.subNet.(networkNames{net}).Load{1,n}) && nnz(Demand(:,n)) == 0
                nDem = Plant.subNet.(networkNames{net}).Load{1,n}; %The load that this node meets demand for
                for j = 1:1:length(nDem)
                    Demand(:,n) = Demand(:,n) + OptData.Demand.E(1:nS,nDem(j));
                end 
            end 
        end 
          
        PwrTrans = zeros(nS,nN);
        losses = zeros(nS,nN);
        checklosses = zeros(nS,nG);
        for i = 1:1:length(Plant.subNet.Electrical.lineNumber)%Gathering the correct line loss associated with the direction of the power on the lines
            pos = strfind(lineNames{i},'_');
            for j = 1:1:length(Plant.subNet.Electrical.nodes)
                if strcmp(lineNames{i}(1:pos(1)-1),Plant.subNet.Electrical.nodes{j}(1))
                    for t = 1:1:nS
                         PwrTrans(t,j) = PwrTrans(t,j) - OptData.GeneratorState(t,i+nG);
                         if OptData.GeneratorState(t,i+nG) < -50
                            losses(t,j) = losses(t,j) + OptData.TransLoss(t,i+length(Plant.subNet.Electrical.lineNumber));
                            checklosses(t,i) = OptData.TransLoss(t,i+length(Plant.subNet.Electrical.lineNumber))/(-OptData.GeneratorState(t,i+nG));
                         end
                    end
                elseif strcmp(lineNames{i}(pos(2)+1:end),Plant.subNet.Electrical.nodes{j}(1))
                    for t = 1:1:nS-1
                        PwrTrans(t,j) = PwrTrans(t,j) + OptData.GeneratorState(t,i+nG);
                        if OptData.GeneratorState(t,i+nG) > 50
                            losses(t,j) = abs(losses(t,j) + OptData.TransLoss(t,i));
                            checklosses(t,i) = OptData.TransLoss(t,i)/(OptData.GeneratorState(t,i+nG));
                        end 
                    end 
                end 
            end 
            plot(checklosses(:,i))
        end 
        figure(1)
        nodalEnergyBalance =  Generated + PwrTrans - losses - Demand; %Demand already subtracted above           
        plot(checklosses(:,:))
        saveas(gcf, fullfile(folder, 'ElossesN(allcomb).jpg'))
        %testing graphic
        
%         %Create Stacked Floating boxes
% %         y = [1.5 2 3; 4.5 5 6; 7.5 8 9; 10.5 11 12];
% %         f = bar(y,'stacked')
% %         f(1).FaceColor = 'w'
% %         f(1).EdgeColor = 'w'       
%         Transfer = OptData.GeneratorState(:,nG+1:nG+nL);
%         linepos = strfind(lineNames{1},'Plant') + length('Plant'); %same position for every line
%         for i = 1:1:nL
%             for j = 1:1:length(linepos)
%                 if j == 1
%                     line(i,j) = str2num(strcat(lineNames{i}(linepos(1,j))));
%                 else
%                     line(i,j) = str2num(strcat(lineNames{i}(linepos(1,j):end)));
%                 end
%             end
%             ymax = i;
%             x = linspace(line(i,1),line(i,2));
%             xlim = [min(x),max(x)];
%             ylim = [i-0.5, i+0.5];
%             Hlims(i,1:3) = [{xlim},{ylim},ymax];
%         end
%         for t = 1:1:nS
%             for i = 1:1:nL
%                 ymax = Hlims{i,3};
%                 ylim = Hlims{i,2};
%                 for H = 1:1:3
%                     xlim = Hlims{i,1}+(nN*(H-1));
%                     if Transfer(t+(H-1),i) > 0
%                         plot(xlim, [1 1]*ymax, '.-g');
%                     elseif Transfer(t+(H-1),i) < 0 
%                         plot(xlim, [1 1]*ymax, '.-r');
%                     elseif Transfer(t+(H-1),i) == 0
%                          plot(xlim, [1 1]*ymax, '.--y')
%                     end
%                     hold on
%                 end
%                 axis([0 nN*3+1   0  nL+1])
%             end   
% %             axis([0 nN*3+1   0  nL+1])
%         end 
%         hold off

        %% Start Electricity Graphing
        for i = 1:1:nN %number of nodes
            plot(Date,nodalEnergyBalance(:,i),'b'); %blue
            hold on
            %Blue = Generated; Green = Power into node; yellow = power out of node; red line = Power out of node; Red stars line = Demand
            if sum(PwrTrans(:,i)) > 0
                area(Date,[Generated(:,i), PwrTrans(:,i), losses(:,i)])
            else
                area(Date,[Generated(:,i), losses(:,i),PwrTrans(:,i)])
            end
            plot(Date,Demand(:,i),'-*r')
            hold off
            saveas(gcf, fullfile(folder, sprintf('EnergyBalanceN%d.jpg', i)))
            
            yyaxis left
            plot(Date,PwrTrans(:,i),'b') %Power Transfers
            hold on
            yyaxis right
            plot(Date,(abs(PwrTrans(:,i))./PwrTrans(:,i)).*losses(:,i),'-.r') %actual
            hold off
            saveas(gcf, fullfile(folder, sprintf('Elossescomp N%d.jpg', i)))
            cla('reset') %resets axis
        end
        plot(Date(2:end),nodalEnergyBalance(2:nS,:)); %First and last point bad; extremely high
        saveas(gcf, fullfile(folder, 'EnergyBalanceN(allcomb).jpg'))
        
        plot(Date,losses)
        saveas(gcf, fullfile(folder, 'Elosses(allcomb).jpg'))
        
        for i = 1:1:nG
            plot(Date,Hist.PowerGen(:,i),'g')
            hold on
            plot(Date(:,1),OptData.GeneratorState(1:nS,i),'--r')
            hold off
            saveas(gcf, fullfile(folder, sprintf('EnergyGenD%d.jpg', i)))
            
            Hper(:,i) = Hist.PowerGen(:,i)./Hist.PowerFlow(:,i);
            Aper(:,i) = OptData.GeneratorState(1:nS,i)./OptData.PowerFlow(1:nS,i);
            plot(Hper(:,i),'r')
            hold on
            plot(Aper(:,i),'b')
            hold off
            saveas(gcf, fullfile(folder, sprintf('EnergyconvdiffD%d.jpg', i)))
            
%             for j = 1:1:length(Hper(:,i))
%                 plot(Hist.PowerGen(j,i),Hper(j,i),'*g');
%                 hold on
%                 plot(OptData.GeneratorState(1:nS,i),Aper(j,i),'or')
%             end
            hold off
            saveas(gcf,fullfile(folder,sprintf('EHisGenHisPF AGenAPF D%d.jpg', i)))
            
        end
        
    elseif strcmp((networkNames{net}),'Hydro') %Hydro mass, SOC, and Outflow balances
       %% For optimization Data
       InFlow = OptData.InFlow(1:nS,:);
       OutFlow = OptData.OutFlow(1:nS,:);
       PowerFlow = OptData.PowerFlow(1:nS,:);
       SpillFlow = OutFlow - PowerFlow;
       storedFlow = (InFlow - OutFlow);
       ChangeSOC = (InFlow - OutFlow)/(dt*12.1);
       HistChangeSOC = (Hist.InFlow-Hist.OutFlow)/(dt*12.1);
  
       %% Start Hydro Graphing
       for i = 1:1:nG
           h = area(Date,[PowerFlow(:,i), SpillFlow(:,i)]); %PowerFlow = purple; spillFlow = yellow;
           h(1).FaceColor = [1,.90,0];
           h(2).FaceColor = [0,1,.70];
           hold on
           plot(Date,InFlow(:,i),'--r'); %red dashed line
           plot(Date,OutFlow(:,i),'b'); %blue
           hold off
           saveas(gcf, fullfile(folder, sprintf('WaterBalanceD%d.jpg', i)))
           
           h = area(Date,[Hist.PowerFlow(:,i), Hist.SpillFlow(:,i)]); %PowerFlow = purple; spillFlow = yellow;
           h(1).FaceColor = [1,.90,0];
           h(2).FaceColor = [0,1,.70];
           hold on
           plot(Date,Hist.InFlow(:,i),'--r'); %red dashed line
           plot(Date,Hist.OutFlow(:,i),'b'); %blue
           hold off
           saveas(gcf, fullfile(folder, sprintf('WaterHistBalanceD%d.jpg', i)))
           
           h = area(Date,PowerFlow(:,i)); 
           h(1).FaceColor = [1,.90,0];
           hold on
           plot(Date,Hist.PowerFlow(:,i),'k');
           hold off
           saveas(gcf, fullfile(folder, sprintf('WaterpowercD%d.jpg', i)))
           
           plot(Date,PowerFlow(:,i).*Aper(:,i),'*r');
           hold on
           plot(Date,OptData.GeneratorState(1:nS,i),'b');
           plot(Date,Hist.PowerFlow(:,i).*Hper(:,i),'*y');
           plot(Date,Hist.PowerGen(:,i),'c');
           hold off
           saveas(gcf, fullfile(folder, sprintf('EpowerFlowcD%d.jpg', i)))
           
           h = plot(Date,SpillFlow(:,i)); 
           h(1).Color = [0 0.4 0];
           hold on
           plot(Date,Hist.SpillFlow(:,i),'g');
           hold off
           saveas(gcf, fullfile(folder, sprintf('WaterspillcD%d.jpg', i)))
           
%            plot(Date,OutFlow(:,i),'b'); %blue
%            hold on
           plot(Date,Hist.OutFlow(:,i),'g'); %blue
%            hold off
           saveas(gcf, fullfile(folder, sprintf('WaterOutFlowsD%d.jpg', i)));
           
%            plot(Date,InFlow(:,i),'k'); %red dashed line
%            hold on
           plot(Date,Hist.InFlow(:,i),'r'); %red dashed line
%            hold off
           saveas(gcf, fullfile(folder, sprintf('WaterInFlowsD%d.jpg', i)));

%            area(Date,ChangeSOC(:,i)); %Blue
%            hold on
           plot(Date,HistChangeSOC(:,i),'g'); %greem
%            hold off
           saveas(gcf, fullfile(folder, sprintf('WaterSOCD%d.jpg', i)));
           
           plot(Date,InFlow(:,i),'--b'); %red dashed line
           hold on
           plot(Date,OutFlow(:,i),'c'); %blue
           plot(Date,Hist.InFlow(:,i),'--g');
           h = plot(Date,Hist.OutFlow(:,i),'g');
           h.Color = [0 .5 0];
           hold off
           saveas(gcf, fullfile(folder, sprintf('WaterHistvsBalanceD%d.jpg', i)));
       end 
       
       plot(Date,ChangeSOC)
       saveas(gcf, fullfile(folder, 'WaterBalanceN(allcomb).jpg')) 
       
       plot(Date,HistChangeSOC)
       saveas(gcf, fullfile(folder, 'WaterBalanceHistorN(allcomb).jpg')) 
       
       for i=1:1:nG
            plot(Date,Hper(:,i)-Aper(:,i),'r')
            hold on
            plot(Date,zeros(nS),'g'); %zero line
            me = mean(Hist.PowerFlow(:,i));
            plot(Date,(Hist.PowerFlow(:,i)-me).*100,'b')
            plot(Date,zeros(nS),'g'); %zero line
%             plot(Date,PowerFlow(:,i)-Hist.PowerFlow(:,i),'b')
           hold off
           saveas(gcf, fullfile(folder, sprintf('WED%d.jpg', i)))
           
           flowOver(:,i) = (OptData.GeneratorState(1:nS,i)-OptData.Demand.E(1:nS,i))*Plant.Generator(i).QPform.Stor.Power2Flow;
           justPF(:,i) = PowerFlow(:,i) - flowOver(:,i);
           for j = 1:1:length(flowOver(:,i))
               if flowOver(j,i) > 0
                   PFplus(j,i) = PowerFlow(j,i) - flowOver(j,i);
                   PFunder(j,i) = 0;
               elseif flowOver(j,i) < 0
                   PFplus(j,i) = 0;
                   PFunder(j,i) = PowerFlow(j,i) - flowOver(j,i);
               end
           end
           plot(flowOver(:,i),'--y')
           hold on
           plot(Hist.PowerFlow(:,i),'--k')
           plot(PowerFlow(:,i),'--c')
           plot(justPF(:,i),'g')
           hold off
           saveas(gcf, fullfile(folder, sprintf('WENoExcessPF D%d.jpg', i)))
           
            plot(flowOver(:,i),'--y')
            hold on
            plot(Hist.PowerFlow(:,i),'--k')
            plot(PFunder(:,i),'r')
            plot(PFplus(:,i),'g')
            hold off
            saveas(gcf, fullfile(folder, sprintf('WENoOverUnderPF D%d.jpg', i)))
       end
           
      for i = 1:1:nG 
           cla('reset')
           left_color = [1 0 0];
           right_color = [0 0 1];
           set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
           yyaxis left
           for j = 1:1:nS
               if Hist.SpillFlow(j,i)>0
                    plot(Date(j,1),SpillFlow(j,i)/Hist.SpillFlow(j,i),'*r')
               end
               hold on
           end
           yyaxis right
           plot(Date(:,1),PowerFlow(:,i)/Hist.PowerFlow(:,i),'b');
           saveas(gcf, fullfile(folder, sprintf('WratioOptvsHistD%d.jpg', i)));
           hold off
           cla('reset') %resets axis
      end 
      
    end 
end
