function [T , Y] = solverFixedStep(RunTime,Step,IC)
n = ceil(RunTime/Step);
T = (0:Step:n*Step)';
Y = zeros(n+1,length(IC));
Y(1,:) = IC';
for t = 1:1:n
    dY = RunBlocks((t-1)*Step,Y(t,:)');
    Y(t+1,:) = Y(t,:) + Step*dY';
end