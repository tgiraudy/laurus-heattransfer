function [ R ] = naturalConvectionVerticalClosedEnclosure( Tinv,Tout,Geometry,MaterialProperties )
%NATURALCONVECTIONVERTICALCLOSEDENCLOSURE Summary of this function goes here
%   Detailed explanation goes here
airProperties = AirProperties();

gbetanu2 = airProperties.gbetanu2;
Pr = airProperties.Pr;
k = airProperties.k;

L = Geometry.L; % altura
b  = Geometry.b; % espaciado entre placas

A = Geometry.A;

Gr = gbetanu2 * abs(Tinv - Tout) * delta^3;
Ra = Gr*Pr;

LdeltaRatio = L/delta;

if LdeltaRatio > 2 && LdeltaRatio < 10 && Pr<10 && Ra < 1e10
    Nu = 0.22 * (LdeltaRatio)^(-1/4) * (Pr/(0.2 + Pr)*Ra)^.28;
elseif LdeltaRatio > 1 && LdeltaRatio < 2 &&  Pr > 1e-3 && LdeltaRatio < 1e5 && Ra * Pr / (0.2+Pr)>1e3
    Nu = 0.18 * (Pr / (0.2 + Pr) * Ra)^.29;
else
    disp('Correlacion para conveccion natural fuera de rango')
end

h = Nu*k/delta;

R = 1/(h*A);


end

