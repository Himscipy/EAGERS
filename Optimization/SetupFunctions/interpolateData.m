function interpData = interpolateData(timestamp,Variability)
%This function will interpolate or average to create data points in-line with the specifid frequency in s
%If the resolution of optimization is faster than the stored data it will interpolate and add noise
%If the resolution is larger, it will average the data within each step of the optimization
%% !!current limitation is that all saved data streams have the same timestamp 
%% !!all saved data has constant sampling frequency
global TestData
Tfreq = (timestamp(2) - timestamp(1))*3600*24;
Xi = nnz(TestData.Timestamp<=timestamp(1));
Xf = nnz(TestData.Timestamp<=timestamp(end));
r = (TestData.Timestamp(2)-TestData.Timestamp(1))*24*3600/Tfreq;% # of points that must be created between datum
lengthTD = length(TestData.Timestamp);
NumSteps = length(timestamp);
%% create noise vector
if r>0 
    z = randn(NumSteps,1);
    N = zeros(NumSteps,1);
    N(1) = 0;
    N(2) = (rand(1,1)-.5)*Variability; %the value .04 makes the final noise signal peaks = variability with this strategy.
    b= sign(N(2)-N(1));
    c = N(2);
    for n = 3:length(N)
        %if the noise is increasing the probability is such that the noise should continue to increase, but as the noise approaches the peak noise magnitude the probability of switching increases.
        if (c>0 && N(n-1)>0) || (c<0 && N(n-1)<0)
            a = 2*(Variability - abs(N(n-1)))/Variability; %constant 2 = 97.7% probability load changes in same direction, 1.5 = 93.3%, 1 = 84.4%, .5 = 69.1%
        else a = 2;
        end
        if abs(z(n))>.68 %only changes value 50% of the time
             c = b*(z(n)+a); % c is positive if abs(Noise)is increasing, negative if decreasing
             b = sign(c);
             N(n) = N(n-1)+c*Variability; 
        else N(n) = N(n-1);
        end
    end
    Noise = N;% scaled noise to a signal with average magnitude of 1
end
interpData.Timestamp = timestamp;
D = datevec(interpData.Timestamp(1));
h_of_y = linspace(8760/length(TestData.Weather.Tdb),8760,length(TestData.Weather.Tdb))';% Hour since start of year
S = {'Tdb';'Twb';'irradDireNorm';};
a = [0;h_of_y];
t = mod(24*(interpData.Timestamp - datenum([D(1),1,1])),8760);% Hour since start of year
for i = 1:1:length(S)
    b = [TestData.Weather.(S{i})(1); TestData.Weather.(S{i})];
    interpData.Weather.(S{i}) = interp1(a,b,t); 
end
F = fieldnames(TestData);
F = F(~strcmp('Timestamp',F));
F = F(~strcmp('Weather',F));
%% take data exactly, or interpolate if necessary
if abs(r-1)<1e-8
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                interpData.(F{j}).(S{i}) = TestData.(F{j}).(S{i})(Xi:Xf,:);
            end
        else
            interpData.(F{j}) = TestData.(F{j})(Xi:Xf);
        end
    end
elseif r<1 %extra datum, average points in between
    x1 = Xi;
    for t = 1:1:NumSteps
        x2 = nnz(TestData.Timestamp<=timestamp(t));
        for j = 1:1:length(F)
            if isstruct(TestData.(F{j}))
                S = fieldnames(TestData.(F{j}));
                for i = 1:1:length(S)
                    for k = 1:1:length(TestData.(F{j}).(S{i})(x1,:))
                        interpData.(F{j}).(S{i})(t,k) = mean(TestData.(F{j}).(S{i})(x1:x2,k));
                    end
                end
            else
                for k = 1:1:length(TestData.(F{j})(x1,:))
                    interpData.(F{j})(t,k) = mean(TestData.(F{j})(x1:x2,k));
                end
            end
        end
        x1=x2+1;
    end
elseif r>1 %interpolate between timesteps
    for j = 1:1:length(F)
        if isstruct(TestData.(F{j}))
            S = fieldnames(TestData.(F{j}));
            for i = 1:1:length(S)
                interpData.(F{j}).(S{i}) = zeros(NumSteps,length(TestData.(F{j}).(S{i})(1,:))); 
            end
            x1 = Xi;
            x2 = x1+1;
            for t = 1:1:NumSteps
                r0 = (mod(t-1,r))/r;
                if TestData.Timestamp(x2)<=interpData.Timestamp(t)
                    x1 = x1+1;
                    x2 = min(x1+1,lengthTD);
                end
                for i = 1:1:length(S)
                    interpData.(F{j}).(S{i})(t,:) = ((1-r0)*TestData.(F{j}).(S{i})(x1,:)+r0*TestData.(F{j}).(S{i})(x2,:))*(1+Noise(t));
                end
            end
        else
            interpData.(F{j}) = zeros(NumSteps,length(TestData.(F{j})(1,:))); 
            x1 = Xi;
            x2 = x1+1;
            for t = 1:1:NumSteps
                r0 = (mod(t-1,r))/r;
                if TestData.Timestamp(x2)<=interpData.Timestamp(t)
                    x1 = x1+1;
                    x2 = min(x1+1,lengthTD);
                end
                interpData.(F{j})(t,:) = ((1-r0)*TestData.(F{j})(x1,:)+r0*TestData.(F{j})(x2,:))*(1+Noise(t));
            end
        end
    end
end