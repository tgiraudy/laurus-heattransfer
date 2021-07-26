function R = conductionHorizClosedEnclosures(Tinv,Tout,Geometry,MaterialProperties)

air = Air();

k = air.k;

L = Geometry.L;
A = Geometry.A;

Nu = 1;

h = Nu*k/L;

q = h * (Tinv - Tout) * A;
R = 1/(h*A);
    