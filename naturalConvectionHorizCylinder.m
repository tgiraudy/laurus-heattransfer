function [R] = naturalConvectionHorizCylinder(Tinv,Tout,Geometry,MaterialProperties)

fluid = Air();

gbetanu2 = fluid.gbetanu2;
Pr = fluid.Pr;
k = fluid.k;

D = Geometry.D;
A = Geometry.A;

Gr = gbetanu2 * abs(Tinv - Tout) * D^3;
Ra = Gr*Pr;

if Ra < 1e12
    Nu = (0.6 + 0.387*Ra^(1/6)/((1+(0.559/Pr)^(9/16))^(8/27)))^2;
else
    disp('Correlacion para conveccion natural en cilindros fuera de rango')
end

if Tinv < Tout
    Nu = -Nu;
end
h = Nu*k/D;

q = h * (Tinv - Tout) * A;
R = 1/(h*A);    