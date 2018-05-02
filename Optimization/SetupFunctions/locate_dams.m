function [gen,subnet] = locate_dams(gen,subnet,dt)
n_g = length(gen);
n_s = length(dt);
upriver = {};
downriver = {};
node_index = (1:1:length(subnet.Hydro.lineNames));
for i = 1:1:length(subnet.Hydro.lineNames)
    name = subnet.Hydro.lineNames{i};
    k = strfind(name,'_');
    upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
    downriver(end+1) = {name(k(2)+1:end)};% node names of the downriver node
end
for i = 1:1:n_g
    if strcmp(gen(i).Type,'Hydro Storage')
        n = gen(i).QPform.Hydro.subnetNode;
        segment_dr = subnet.Hydro.lineNumber(strcmp(subnet.Hydro.nodes{n}(1),upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
        subnet.Hydro.DownRiverSegment(n) = segment_dr;
        if any(strcmp(subnet.Hydro.nodes{n}(1),downriver))
            subnet.Hydro.UpRiverNodes(n) = {node_index(strcmp(subnet.Hydro.nodes{n}(1),downriver))};%upstream nodes
            segment_ur = subnet.Hydro.lineNumber(subnet.Hydro.UpRiverNodes{n});
        else
            segment_ur = []; %no upstream nodes
        end
        for t = 1:1:n_s
            %river segments flowing into this node
            for j = 1:1:length(segment_ur)
                segment_time = subnet.Hydro.lineTime(subnet.Hydro.UpRiverNodes{n}(j));
                tt = sum(dt(1:t));
                if tt<segment_time                            
                    %Do nothing; the inflow rate will be updated in update_matrices
                elseif tt>=segment_time && tt<=segment_time+dt(1)%between initial condition & first step
                    frac = (tt-segment_time)/dt(1);%portion of flow from step 1, remainder from  SourceSink + (1-frac)*Inflow(t=1) : subtracted from beq in update matrices
                    subnet.Hydro.frac(subnet.Hydro.UpRiverNodes{n}(j)) = frac;
                end 
            end 
            %water flow out of the node 
            gen(i).QPform.DownRiverSegment = segment_dr;
        end 
    else
        %% add water district here
    end
end