function [J, R] = residualCalculation(T, prevT, bc, elements, timeStep)

    L = 6;
    D = 4;
    V = L*pi*D^2/4/2;

    bAislacion = 0.02;
    bSuelo = 0.2;

    nEle = size(elements, 1);
    nNodEle = size(elements, 2);
    nNod = max(max(elements));
    nDofNod = 1;
    nDofEle = nDofNod*nNodEle;
    nDofTot = nDofNod*nNod;
    dofs = reshape(1:nDofTot, nDofNod, [])';

    Tinf = bc(1);
    T0 = bc(2);
    Tsky = bc(3);
    sunRadiationHeat = bc(4);
    airMass = bc(5);

    solarAbsorptivity = 0.5;
    floorEmissivity = 0.5;
    air = Air();

    K = zeros(nDofTot);
    C = zeros(nDofTot);
    q = zeros(nDofTot, 1);

    q(4) = q(4) + sunRadiationHeat;

    C(1,1) = airMass*air.Cp;

    R1 = 1e-4/(T(1) - T(2)); % Esta resistencai es recontra falopa

    Kc1 = [1/R1 -1/R1
        -1/R1  1/R1];
    K(dofs(elements(1,:), :), dofs(elements(1,:), :)) = K(dofs(elements(1,:), :), dofs(elements(1,:), :)) + Kc1;
    q12 = 1/R1*(T(2) - T(1));
    q(1) = q(1) + q12;

    Geometry.L = L;
    Geometry.Di = D;
    Geometry.b = bAislacion;
    Geometry.theta = pi; %semicylinder
    R2 = naturalConvectionCylinderClosedEnclosures(T(2), T(3), Geometry, []);


    Geometry.D = D;
    Geometry.A = pi*D*L/2;
    Rinf = naturalConvectionHorizCylinder(T(3), T(4), Geometry, [])/2;


    Kc2 = [1/R2 -1/R2
        -1/R2  1/R2 + 1/Rinf];

    qcinf = 1/Rinf*(Tinf - T(3));
    qc23 = 1/R2*(T(3) - T(2));

    q(2) = q(2) + qc23 - q12;
    q(3) = q(3) + qcinf - qc23;

    K(dofs(elements(2,:), :), dofs(elements(2,:), :)) = K(dofs(elements(2,:), :), dofs(elements(2,:), :)) + Kc2;

    Geometry.D = D;
    Geometry.L = L;
    Geometry.A = L*D;
    Geometry.P = 2*L+2*D;
    R3 = naturalConvectionHorizFlatPlate(T(1),T(4),Geometry,[]);


    Geometry.b = bSuelo;
    Geometry.A = L*D;
    R0 = naturalConvectionHorizClosedEnclosures(T(4),T0,Geometry,[]);

    qc0 = 1/R0*(T0 - T(4));
    q(4) = q(4) + qc0;

    Kc3 = [1/R3 -1/R3
        -1/R3  1/R3 + 1/R0];

    Geometry.A = L*D;
    MaterialProperties.emissivity = floorEmissivity;
    hrad = reradiation(T(4), Tsky, Geometry, MaterialProperties);

    Kr3 = [0 0
        0 4*hrad*T(4)^3];

    qr = hrad*(Tsky^2 + T(4)^2)*(Tsky + T(4));
    q(4) = q(4) + qr;

    q = q - C*(T - prevT)/(timeStep*3600);
    
    K(dofs(elements(3,:), :), dofs(elements(3,:), :)) = K(dofs(elements(3,:), :), dofs(elements(3,:), :)) + Kc3 + Kr3;

    J = K + C/(timeStep*3600);
    R = q;
end