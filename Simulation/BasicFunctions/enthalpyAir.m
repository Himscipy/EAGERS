function h = enthalpyAir(Air)
% calculates saturated specific enthalpy (kJ/kg dry air) of a flow of air (HVAC)
w = Air.H2O*18./(MassFlow(Air)-Air.H2O*18);%kg H20/kg dry air
h = (1.006*(Air.T-273.15) + w.*(2501 + 1.86*(Air.T-273.15))); %ASHRAE 2013 fundamentals eq. 32
end%Ends function enthalpyAir