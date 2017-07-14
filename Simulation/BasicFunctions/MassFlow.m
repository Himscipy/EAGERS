function MF = MassFlow(Flow)
spec = fieldnames(Flow);
MF = 0;
for i = 1:1:length(spec)
    switch spec{i}
        case 'CH4'
            MF = MF + Flow.(spec{i})*16;
        case 'CO'
            MF = MF + Flow.(spec{i})*28;
        case 'CO2'
            MF = MF + Flow.(spec{i})*44;
        case 'H2'
            MF = MF + Flow.(spec{i})*2;
        case 'H2O'
            MF = MF + Flow.(spec{i})*18;
        case 'N2'
            MF = MF + Flow.(spec{i})*28;
        case 'O2'
            MF = MF + Flow.(spec{i})*32;
        case 'C'
            MF = MF + Flow.(spec{i})*12;
        case 'NO'
            MF = MF + Flow.(spec{i})*30;
        case 'AR'
            MF = MF + Flow.(spec{i})*40;
        case 'OH'
            MF = MF + Flow.(spec{i})*17;
        case 'H'
            MF = MF + Flow.(spec{i})*1;
        case 'C2H6'
            MF = MF + Flow.(spec{i})*30;
        case 'C3H8'
            MF = MF + Flow.(spec{i})*44;
        case 'C4H10'
            MF = MF + Flow.(spec{i})*58;
        case 'C5H12'
            MF = MF + Flow.(spec{i})*72;
        case 'C6H14'
            MF = MF + Flow.(spec{i})*86;
        case 'CH3OH'
            MF = MF + Flow.(spec{i})*32;
        case 'C2H5OH'
            MF = MF + Flow.(spec{i})*46;
        case 'C6H6'
            MF = MF + Flow.(spec{i})*78;
    end
end
end%Ends function MassFlow