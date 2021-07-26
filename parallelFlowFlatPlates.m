function [q,R] = parallelFlowFlatPlates(Tinv,Tout,Geometry,MaterialProperties)

air = Air();

k = air.k;
Pr = air.Pr;
U = air.U
rho = air.rho;
mu = air.mu;

L = Geometry.L;
A = Geometry.A;

Re = rho * U * L / mu;

if Re < 5e5
    Nu = 0.664*Re^0.5*Pr^(1/3);
elseif Re < 1e7 && Pr > 0.6 && Pr < 60
    Nu = 0.037*Re^0.8*Pr^(1/3);
else
    disp('correlacion para placa plana no valida')
end

h = Nu*k/L;

q = h * (Tinv - Tout) * A;
R = 1/(h*A);
