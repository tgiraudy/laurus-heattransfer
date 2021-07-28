function [R] = naturalConvectionCylinderClosedEnclosures(Tinv, Tout, Geometry, MaterialProperties)


airProperties = AirProperties();

gbetanu2 = airProperties.gbetanu2;
Pr = airProperties.Pr;
k = airProperties.k;

theta = Geometry.theta;
Di = Geometry.Di;
L = Geometry.L;
b = Geometry.b;

characteristicLength = b;

Do = Di + 2*b;

Gr = gbetanu2 * abs(Tinv - Tout) * characteristicLength^3;
Ra = Gr*Pr;

%p 321 kreith heat transfer

adimCheck = (log(Do./Di)./(b.^.75.*(Do.^(-3/5) + Di.^(-3/5)).^(5/4))).^4.*Ra;
% if Pr >= .7 && Pr <6000 && adimCheck < 1e7 && adimCheck > 1e1;
    kEff = k*0.386*(log(Do./Di)./(b.^.75.*(Do.^(-3/5) + Di.^(-3/5)).^(5/4))).*(Ra.^(1/4))*((Pr/(0.861+Pr)).^(1/4));
% else
%     disp('Correlacion para conveccion natural fuera de rango')
% end


if Tinv < Tout
    kEff = -kEff;
end


% resistencia medio cilindro
R = log(Do./Di)./(theta*kEff*L);
% calor
end