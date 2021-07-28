function [J, R] = residualCalculation(Titer, bc)
        
        bc(1) = Tinf;
        bc(2) = T0;
        bc(3) = Tsky;
        bc(4) = sunRadiationHeat;
        K = zeros(nDofTot);
        C = zeros(nDofTot);
        q = zeros(nDofTot, 1);
        
        q(4) = q(4) + sunRadiationHeat;
        
        C(1,1) = airMass*air.Cp;
        
        R1 = 0.1;
        
        Kc1 = [1/R1 -1/R1
              -1/R1  1/R1];
        K(dofs(elements(1,:), :), dofs(elements(1,:), :)) = K(dofs(elements(1,:), :), dofs(elements(1,:), :)) + Kc1;   
        q12 = 1/R1*(Titer(2) - Titer(1));
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
          
        qcinf = 1/Rinf*(Tinf - Titer(3));
        qc23 = 1/R2*(Titer(3) - Titer(2)); 
        
        q(2) = q(2) + qc23 - q12;
        q(3) = q(3) + qcinf - qc23; 
        
        K(dofs(elements(2,:), :), dofs(elements(2,:), :)) = K(dofs(elements(2,:), :), dofs(elements(2,:), :)) + Kc2;
        
        Geometry.D = D;
        Geometry.L = L;
        Geometry.A = L*D;
        Geometry.P = 2*L+2*D;
        R3 = naturalConvectionHorizFlatPlate(Titer(1),Titer(5),Geometry,[]);
        
        
        Geometry.b = bSuelo;
        Geometry.A = L*D;
        R0 = naturalConvectionHorizClosedEnclosures(Titer(5),Titer(6),Geometry,[]);
        
        qc0 = 1/R0*(T0 - T(4));
        q(4) = q(4) + qc0;
               
        Kc3 = [1/R3 -1/R3
              -1/R3  1/R3 + 1/R0];
        
        Geometry.A = L*D;
        MaterialProperties.emissivity = floorEmissivity;
        hrad = reradiation(Titer(5), Tsky, Geometry, MaterialProperties);
            
        Kr3 = [0 0
               0 4*hrad*Titer(4)^3];
        
        qr = hrad*(Tsky^2 + Titer(4)^2)*(Tsky + Titer(4));
        q(4) = q(4) + qr;
        
        K(dofs(elements(3,:), :), dofs(elements(3,:), :)) = K(dofs(elements(3,:), :), dofs(elements(3,:), :)) + Kc3 + Kr3;
        
        J = K + C/timeStep;
        R = q;
end