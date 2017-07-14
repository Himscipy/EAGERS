function NewAir = MixAir(Air1,Air2)
%combine 2 air streams and find resulting temperature that mathces enthalpy
Air1_mass = MassFlow(Air1) - Air1.H2O*18; %dry air mass flow
Air2_mass = MassFlow(Air2) - Air2.H2O*18; %dry air mass flow
H1 = enthalpyAir(Air1)*Air1_mass;
H2 = enthalpyAir(Air2)*Air2_mass;
NewAir.H2O = Air1.H2O + Air2.H2O;
NewAir.N2 = Air1.N2 + Air2.N2;
NewAir.O2 = Air1.O2 + Air2.O2;
NewAir.AR = Air1.AR + Air2.AR;
NewAir.CO2 = Air1.CO2 + Air2.CO2;
NewAir.T = (Air1.T*Air1_mass + Air2.T*Air2_mass)/(Air1_mass+Air2_mass);
error = 1;
while error>1e-3
    error = (enthalpyAir(NewAir)*(Air1_mass + Air2_mass) - H1 - H2);
    NewAir.T = NewAir.T - error/(1.5*(Air1_mass+Air2_mass));
end
end%Ends function MixAir