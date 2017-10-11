%% create graphs for comparison from Matlab
global Plant InOut 
xi = []; %Start date
xf = []; %End date
nG = length(Plant.Generator); %number of Generators
OptData = Plant.Dispatch;
networkNames = fieldnames(Plant.subNet);

%Calculate balances for graphing
for net = 1:1:length(networkNames)
    if strcmp((networkNames{net}),'Electrical')
        nodal_sumGen = zeros(1,length(Plant.subNet.(networkNames{net}).nodes));
        z = 1;
        for i = 1:1:nG
            n = Plant.Generator(i).QPform.(networkNames{net}).subnetNode; %Electric node number
            nodal_sumGen(:,n) = nodal_sumGen(:,n) + OptData.GeneratorState(:,i); %Sum all power Generated at that node
            if ~isempty(Plant.Network(i).(networkNames{net}).Load)
                Demand(n) = n;
                z = z+1;
            end
        end

        for i = 1:1:length(Plant.subNet.Electrical.nodes)
            for t = 1:1:length(Plant.Dispatch.Trans(:,10))
                if Plant.Dispatch.GeneratorState(t,nG+i) > 0
                    losses(:,i) = Plant.Trans(t,i);
                elseif Plant.Dispatch.GeneratorState(t,nG+i) > 0
                    losses(t,i) = Plant.Trans(t,i+length(Plant.subNet.Electrical.nodes));
                else
                    losses(t,i) = 0;
                end
            end

            if i == Plant.subNet.line
                nodalEnergyBalance(t,i) = nodalEnergyBalance(t,i) + nodal_sumGen(t,i) + Plant.Dispatch.GeneratorState(t,ln(j)) - losses(t,i);
                totalEnergy(t,1)

            end
        end 

        for i = 1:1:length(Demand) %subtract Demand from node
            node = Demand(i);
            nodalEnergyBalance(:,node) = nodalEnergyBalance(:,node) - Plant.Dispatch.Demand.E(t,node);
        end
    elseif strcmp((networkNames{net}),'Hydro')
        
       %Outflow 
       %Mass balance
       %Storage Balance
        
    end
    


    %Start Graphing
    for 
        
        
    end
    
end
