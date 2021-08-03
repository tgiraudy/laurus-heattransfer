function R = naturalConvectionVerticalFlatPlate(Tinv,Tout,Geometry,MaterialProperties)

airProperties = AirProperties();

gbetanu2 = airProperties.gbetanu2;
Pr = airProperties.Pr;
k = airProperties.k;
D = Geometry.D;
L = Geometry.L;
A = Geometry.A;
P = Geometry.P;

characteristicLength = A/P;

Gr = gbetanu2 * abs(Tinv - Tout) * characteristicLength^3;
Ra = Gr*Pr;


if Ra > 10 && Ra < 1e8
    Nu = 0.68*Pr^.5*(Gr^.25)/(0.952+Pr)^.25;
elseif Gr>1e9
    Nu = 0.13.*Ra.^(1/3);
else
    Nu = 0.13.*Ra.^(1/3);
    disp('Correlacion para conveccion natural fuera de rango')
end

if Tinv < Tout
    Nu = -Nu;
end

h = Nu*k/characteristicLength;

R = 1/(h*A);
    