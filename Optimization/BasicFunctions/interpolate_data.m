function interp_data = interpolate_data(test_data,freq,variability)
%This function will interpolate or average to create data points in-line with the specifid frequency (Tfreq) in s
%It will add random variability to any interpolated data points
%If the resolution of optimization is faster than the stored data it will interpolate and add noise
%If the resolution is larger, it will average the data within each step of the optimization
%% !!current limitation is that all saved data streams have the same timestamp 
%% !!all saved data has constant sampling frequency
interp_data = [];
if isfield(test_data,'RealTimeData')
    freq2 = (test_data.Timestamp(2) - test_data.Timestamp(1))*24*3600;
    if abs(freq - freq2)<1e-4 %dont re-build historian
        interp_data = test_data.RealTimeData;
    else
        test_data = rmfield(test_data,'RealTimeData');
        if isfield(test_data,'HistProf')
            test_data = rmfield(test_data,'HistProf');
        end
    end
end
if isempty(interp_data)
    timestamp = (test_data.Timestamp(1):freq/24/3600:test_data.Timestamp(end))';
    r = (test_data.Timestamp(2)-test_data.Timestamp(1))*24*3600/freq;% # of points that must be created between datum
    length_td = length(test_data.Timestamp);
    num_steps = length(timestamp);
    %% create noise vector
    if r>0 
        z = randn(num_steps,1);
        noise = zeros(num_steps,1);% scaled noise to a signal with average magnitude of 1
        noise(1) = 0;
        noise(2) = (rand(1,1)-.5)*variability; %the value .04 makes the final noise signal peaks = variability with this strategy.
        b= sign(noise(2)-noise(1));
        c = noise(2);
        for n = 3:length(noise)
            %if the noise is increasing the probability is such that the noise should continue to increase, but as the noise approaches the peak noise magnitude the probability of switching increases.
            if (c>0 && noise(n-1)>0) || (c<0 && noise(n-1)<0)
                a = 2*(variability - abs(noise(n-1)))/variability; %constant 2 = 97.7% probability load changes in same direction, 1.5 = 93.3%, 1 = 84.4%, .5 = 69.1%
            else
                a = 2;
            end
            if abs(z(n))>.68 %only changes value 50% of the time
                 c = b*(z(n)+a); % c is positive if abs(Noise)is increasing, negative if decreasing
                 b = sign(c);
                 noise(n) = noise(n-1)+c*variability; 
            else
                noise(n) = noise(n-1);
            end
        end
    end
    interp_data.Timestamp = timestamp;
    f = fieldnames(test_data);
    f = f(~strcmp('Timestamp',f));
    f = f(~strcmp('HistProf',f));
    %% take data exactly, or interpolate if necessary
    if abs(r-1)<1e-8
        for j = 1:1:length(f)
            if isstruct(test_data.(f{j}))
                s = fieldnames(test_data.(f{j}));
                for i = 1:1:length(s)
                    if isnumeric(test_data.(f{j}).(s{i}))
                        interp_data.(f{j}).(s{i}) = test_data.(f{j}).(s{i});
                    end
                end
            else
                if isnumeric(test_data.(f{j}))
                    interp_data.(f{j}) = test_data.(f{j});
                end
            end
        end
    elseif r<1 %extra datum, average points in between
        x1 = 1;
        for t = 1:1:num_steps
            x2 = nnz(test_data.Timestamp<=timestamp(t));
            for j = 1:1:length(f)
                if isstruct(test_data.(f{j}))
                    s = fieldnames(test_data.(f{j}));
                    for i = 1:1:length(s)
                        if isnumeric(test_data.(f{j}).(s{i}))
                            for k = 1:1:length(test_data.(f{j}).(s{i})(x1,:))
                                interp_data.(f{j}).(s{i})(t,k) = mean(test_data.(f{j}).(s{i})(x1:x2,k));
                            end
                        end
                    end
                elseif isnumeric(test_data.(f{j}))
                    for k = 1:1:length(test_data.(f{j})(x1,:))
                        interp_data.(f{j})(t,k) = mean(test_data.(f{j})(x1:x2,k));
                    end
                end
            end
            x1=x2+1;
        end
    elseif r>1 %interpolate between timesteps
        for j = 1:1:length(f)
            if isstruct(test_data.(f{j}))
                s = fieldnames(test_data.(f{j}));
                for i = 1:1:length(s)
                    if isnumeric(test_data.(f{j}).(s{i}))
                        interp_data.(f{j}).(s{i}) = zeros(num_steps,length(test_data.(f{j}).(s{i})(1,:))); 
                    end
                end
                x1 = 1;
                x2 = x1+1;
                for t = 1:1:num_steps
                    r0 = (mod(t-1,r))/r;
                    if test_data.Timestamp(x2)<=interp_data.Timestamp(t)
                        x1 = x1+1;
                        x2 = min(x1+1,length_td);
                    end
                    for i = 1:1:length(s)
                        if isnumeric(test_data.(f{j}).(s{i}))
                            interp_data.(f{j}).(s{i})(t,:) = ((1-r0)*test_data.(f{j}).(s{i})(x1,:)+r0*test_data.(f{j}).(s{i})(x2,:))*(1+noise(t));
                        end
                    end
                end
            else
                if isnumeric(test_data.(f{j}))
                    interp_data.(f{j}) = zeros(num_steps,length(test_data.(f{j})(1,:))); 
                end
                x1 = 1;
                x2 = x1+1;
                for t = 1:1:num_steps
                    r0 = (mod(t-1,r))/r;
                    if test_data.Timestamp(x2)<=interp_data.Timestamp(t)
                        x1 = x1+1;
                        x2 = min(x1+1,length_td);
                    end
                    if isnumeric(test_data.(f{j}))
                        interp_data.(f{j})(t,:) = ((1-r0)*test_data.(f{j})(x1,:)+r0*test_data.(f{j})(x2,:))*(1+noise(t));
                    end
                end
            end
        end
    end
end
end%Ends function interpolate_data