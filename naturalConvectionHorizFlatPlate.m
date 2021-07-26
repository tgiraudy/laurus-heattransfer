function R = naturalConvectionHorizFlatPlate(Tinv,Tout,Geometry,MaterialProperties)

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


if Ra > 1e5 && Ra < 1e11
    Nu = 0.27.*Ra.^(1/4);
else
    disp('Correlacion para conveccion natural fuera de rango')
end

h = Nu*k/characteristicLength;

R = 1/(h*A);
    