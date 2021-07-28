function [hr] = reradiation(Tfloor, Tout, Geometry, MaterialProperties)

%%%%%%
% inputs
% TFloor: Temperatura del suelo
% TOut: Temperatura del ambiente al que radia (temperatura exterior)
% Geometry: para este caso solo area del suelo (.A);
% MaterialProperties solo necesita emisividad (.emissivity)
% 
% output: resistencia
%%%%%%
c2k = @(c) c + 273.15;

emissivity = MaterialProperties.emissivity;
surfaceArea = Geometry.A;

sb = 5.67e-8; % stefan - boltzmann (W/m^2K^4)

Tfloor = c2k(Tfloor); %K
Tout = c2k(Tout); %K

% q = surfaceArea*sb*emissivity*(Tfloor.^4 - Tout.^4);

hr = surfaceArea*sb*emissivity;

end