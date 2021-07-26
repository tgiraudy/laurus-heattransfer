function [R] = ventilationLoss(Tinv, Tout, Geometry, MaterialProperties)

air = Air();

airFlow = Geometry.airFlow;
V = Geometry.V;

massFlow = air.rho*airFlow;

qLost = massFlow*air.Cp*(Tinv - Tout);

R = (Tinv - Tout)/qLost;

end
