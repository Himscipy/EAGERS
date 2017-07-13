function [LHV, HHV] = HeatingValue(Flow)
spec = fieldnames(Flow);
LHV = 0;
H2O = 0;
for i = 1:1:length(spec)
    switch spec{i}
        case 'CH4'
            LHV = LHV + Flow.(spec{i})*802952.15; %kmol * kJ/kmol
            H2O = H2O + 2*Flow.(spec{i}); %kmol of H2O produced
        case 'CO'
            LHV = LHV + Flow.(spec{i})*305200; %kmol * kJ/kmol
            %no H2O
        case 'CO2'
            %no heating value
            %no H2O
        case 'H2'
            LHV = LHV + Flow.(spec{i})*240424; %kmol * kJ/kmol
            H2O = H2O + Flow.(spec{i});
        case 'H2O'
            %no heating value
            %no H2O
        case 'N2'
            %no heating value
            %no H2O
        case 'O2'
            %no heating value
            %no H2O
        case 'C'
            LHV = LHV + Flow.(spec{i})*393500; %kmol * kJ/kmol
            %no H2O
        case 'NO'
            %no heating value
            %no H2O
        case 'AR'
            %no heating value
            %no H2O
        case 'OH'
            %no heating value
            %no H2O
        case 'H'
            %no heating value
            %no H2O
        case 'C2H6'%ethane
            LHV = LHV + Flow.(spec{i})*1428700; %kmol * kJ/kmol
            H2O = H2O + 3*Flow.(spec{i}); %kmol of H2O produced
        case 'C3H8'%propane
            LHV = LHV + Flow.(spec{i})*2039400; %kmol * kJ/kmol
            H2O = H2O + 4*Flow.(spec{i}); %kmol of H2O produced
        case 'C4H10'%butane
            LHV = LHV + Flow.(spec{i})*2653500; %kmol * kJ/kmol
            H2O = H2O + 5*Flow.(spec{i}); %kmol of H2O produced
        case 'C5H12'%pentane
            LHV = LHV + Flow.(spec{i})*3265200; %kmol * kJ/kmol
            H2O = H2O + 6*Flow.(spec{i}); %kmol of H2O produced
        case 'C6H14'%hexane
            LHV = LHV + Flow.(spec{i})*3875380; %kmol * kJ/kmol
            H2O = H2O + 7*Flow.(spec{i}); %kmol of H2O produced
        case 'CH3OH'%methanol
            LHV = LHV + Flow.(spec{i})*726000; %kmol * kJ/kmol
            H2O = H2O + 2*Flow.(spec{i}); %kmol of H2O produced
        case 'C2H5OH'%ethanol
            LHV = LHV + Flow.(spec{i})*1300000; %kmol * kJ/kmol
            H2O = H2O + 3*Flow.(spec{i}); %kmol of H2O produced
        case 'C6H6'%Benzene
            LHV = LHV + Flow.(spec{i})*3270000; %kmol * kJ/kmol
            H2O = H2O + 3*Flow.(spec{i}); %kmol of H2O produced
    end
end
HHV = LHV + H2O*40660; %latent heat of water = 40660 kJ/kmol
LHV = LHV/NetFlow(Flow); %calculate heating value in kJ/kmol
HHV = HHV/NetFlow(Flow); %calculate heating value in kJ/kmol
end%Ends function HeatingValue