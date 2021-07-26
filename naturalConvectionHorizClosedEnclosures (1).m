function [q, R] = naturalConvectionHorizClosedEnclosures(Tinv,Tout,Geometry,MaterialProperties)

air = Air();

gbetanu2 = air.gbetanu2;
Pr = air.Pr;
k = air.k;

L = Geometry.L;
A = Geometry.A;

Gr = gbetanu2 * abs(Tinv - Tout) * L^3;
Ra = Gr*Pr;

if Ra > 1e4 && Ra < 4e5
    Nu = 0.195*Ra^(1/4);
elseif Ra < 1e8
    Nu = 1 + 1.44*abs(1 - (1708/Ra)) + abs((Ra^(1/3)/18) - 1);
else
    disp('Correlacion para conveccion natural fuera de rango')
end

h = Nu*k/L;

q = h * (Tinv - Tout) * A;
R = 1/(h*A);    