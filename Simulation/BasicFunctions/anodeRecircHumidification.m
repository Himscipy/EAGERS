function r = anodeRecircHumidification(Spec,Flow,WGSeffective,Steam2Carbon,H2Oproduced,r)
CH4 = Spec.CH4*Flow;
dr = 1e-6;
error = 1;
% Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
while abs(error)>1e-6
    COin = Spec.CO*Flow/(1-r) + (Spec.CH4 - WGSeffective*(Spec.CH4+Spec.CO))*Flow*r/(1-r);
    Fuel_H2O = Spec.H2O*Flow/(1-r) + (H2Oproduced - (Spec.CH4 + (Spec.CH4 + Spec.CO)*WGSeffective)*Flow)*r/(1-r);
    S2C = Fuel_H2O/(CH4 + 0.5*COin);
    error = Steam2Carbon - S2C;
    r2 = r+dr;
    COin2 = Spec.CO*Flow/(1-r2) + (Spec.CH4 - WGSeffective*(Spec.CH4+Spec.CO))*Flow*r2/(1-r2);
    H2Oin2 = Spec.H2O*Flow/(1-r2) + (H2Oproduced - (Spec.CH4 + (Spec.CH4 + Spec.CO)*WGSeffective)*Flow)*r2/(1-r2);
    S2C2 = H2Oin2/(CH4 + 0.5*COin2);
    dSdr = (S2C2 - S2C)/dr;
    r = r + error/dSdr;
end
end%Ends function anodeRecircHumidification