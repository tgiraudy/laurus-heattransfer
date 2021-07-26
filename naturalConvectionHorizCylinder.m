function [q, R] = naturalConvectionHorizCylinder(Tinv,Tout,Geometry,MaterialProperties)

fluid = Air();

gbetanu2 = fluid.gbetanu2;
Pr = fluid.Pr;
k = fluid.k;

L = Geometry.L;
A = Geometry.A;

Gr = gbetanu2 * abs(Tinv - Tout) * L^3;
Ra = Gr*Pr;

if Gr > 1e3 && Gr < 1e9
    Nu = 0.53 * (Ra) ^ 0.25 ;
else
    disp('Correlacion para conveccion natural en cilindros fuera de rango')
end

h = Nu*k/L;

q = h * (Tinv - Tout) * A;
R = 1/(h*A);    