%% Ubicacion

Options.Ubicacion = [-34.4 -58.6];
Options.GMT = -3;
Options.Albedo = 1;

%% elements

elements = [1 2   % Conveccion natural interna
            2 3   % Conveccion en enclosure
            3 4   % Conveccion natural externa
            1 5   % Conveccion natural con piso
            5 6   % Conveccion natural piso-suelo
            5 7]; % Reradiacion piso-cielo

nEle = size(elements, 1);
nNodEle = size(elements, 2);
nNod = max(max(elements));
nDofNod = 1;
nDofEle = nDofNod*nNodEle;


Tseed = 25;
Kele = zeros(nDofEle, nDofEle, nEle);
T = Tseed*ones(nNod, 1);
q = zeros(nNod, 1);
        
%% Geometria        
L = 15;
D = 10;
bAislacion = 0.02;
bSuelo = 0.2;
        
k2c = @(T) T - 273.15;        
% si se agregan costados tiene que entrar en el calculo de elementos 1, 2 y
% 3;

%% Propiedades

solarAbsorptivity = 0.5;
floorEmissivity = 0.5;

%% DATA

% separar el vector tiempo en dias

nDays = 3;
% for iDays = 1:nDays
% hour = timeVector;
% day = 196;
% radiationVector = zeros(1, length(hour));
% 
% for i = 1:size(radiationVector, 2)
%     radiationVector(i) = SunRadiation([0 1 0], hour(i), day, Options);
% end
% end

endTime = 72; %h
timeStep = 0.5; % h

for time = 0:timeStep:endTime; %horas


fixed = [4 6 7];

Tinf = TinfData(tIndex);
hr = hrData(tIndex);
Tdp = Tinf; % Temperatura dew point (EN REALIDAD TIENE QUE SER FUNCION DE Tinf y hr);
epsilonSky = 0.711 + 0.56*(Tdp/100) +  0.73*(Tdp^2/100)^2; % berdahl and martin 1984 REVISAR
Tsky = (epsilonSky)^.25*Tinf;
T(fixed) = [Tinf T0 Tsky];
q(5) = radiationVector(tIndex)*L*D*solarAbsorptivity;


%% Resistance calculation

R1 = 1;

Geometry.L = L;
Geometry.Di = D;
Geometry.b = bAislacion;
Geometry.theta = pi; %semicylinder
R2 = naturalConvectionCylinderClosedEnclosures(T(2), T(3), Geometry, []);

Geometry.D = D;
Geometry.A = pi*D*L/2;
R3 = naturalConvectionHorizCylinder(T(3), T(4), Geometry, [])/2;

Geometry.D = D;
Geometry.L = L;
Geometry.A = L*D;
Geometry.P = 2*L+2*D;
R4 = naturalConvectionHorizFlatPlate(T(1),T(5),Geometry,[]);

Geometry.b = bSuelo;
Geometry.A = L*D;
R5 = naturalConvectionHorizClosedEnclosures(T(5),T(6),Geometry,[]); 

Geometry.A = L*D;
MaterialProperties.emissivity = floorEmissivity;
R6 = reradiation(T(5), Tsky, Gometry, MaterialProperties);


R = [R1 R2 R3 R4 R5 R6];

%% Ensamble de matriz

for iEle = 1:nEle
    eleDofs = dofs(elements(iEle,:),:);
    Kele(:,:, iEle) = [1 -1; -1 1]/R(iEle);
    K(eleDofs, eleDofs) = K(eleDofs, eleDofs) + Kele(:,:, iEle);
end


Tc = T(fixed);
qc = q(~fixed);

Kxx = K(~fixed, ~fixed);
Kxc = K(~fixed, fixed);
Kcx = K(fixed, ~fixed);
Kcc = K(fixed, fixed);

Tx = Kxx\(qc - Kxc*Tc);
qx = Kcc*Tc + Kcx*Tx;

q(fixed) = qx;
T(~fixed) = Tx;



end