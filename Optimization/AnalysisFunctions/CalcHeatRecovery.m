function HeatRecovery = CalcHeatRecovery(Dispatch)
global Plant
nG = length(Plant.Generator);
[nS,~] = size(Dispatch);
HeatRecovery = zeros(nS,1);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'CHP Generator') || (strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(Plant.Generator(i).Source,'Heat'))
        Gen = Plant.Generator(i).QPform;
        states = Gen.states(1:nnz(~cellfun('isempty',Gen.states(:,end))),end);
        for t = 1:1:nS
            pow = 0;
            if Dispatch(t,i)>0 && isfield(Gen,'constDemand') && isfield(Gen.constDemand,'H')
                HeatRecovery(t) = HeatRecovery(t) -Gen.constDemand.H;
            end
            Hratio = Gen.output.H(:,end);
            for j = 1:1:length(states)
                if length(Hratio)>1
                    hr = Hratio(j);
                else hr = Hratio;
                end
                HeatRecovery(t) = HeatRecovery(t) + min(Dispatch(t,i)-pow,Gen.(states{j}).ub(end))*hr;
                pow = pow+min(Dispatch(t,i)-pow,Gen.(states{j}).ub(end));
            end
        end
    end
end
HeatRecovery = HeatRecovery/293.1;
end%ends function CalcHeatRecovery