

data = zeros(24,2);

data(:,1) = 0:23;
data(:,2) = [28 28 26 24 23 21 20 20 19 20 21 25 25 25 27 29 29 31 29 28 28  28 28 27];

L = 15;
D = 10;

alfa = 0.3; % albedo de la superficie del invernadero 

% Options.Ubicacion = [10.66, -73.93];
% Options.GMT = -4;
% Options.Albedo = 1.0;
% hour = 0:23;
% day = 40;



Options.Ubicacion = [-34.4 -58.6];
Options.GMT = -3;
Options.Albedo = 1;
hour = 12;
day = 196;
radiationVector = zeros(1, length(hour));

for i = 1:size(radiationVector, 2)
    radiationVector(i) = SunRadiation([0 1 0], hour(i), day, Options);
end

f2c = @(f)(f-32)*5/9;
time = data(:,1);
T2 = data(:,2);

Geometry.A = L*D;

qRadiation = radiationVector*Geometry.A*alfa;


disp(['Incoming Radiative Heat : ' num2str(qRadiation/1000) ' kW'])


Tinv = 21;
Tout = 10;

Geometry.L = 0.02;
Geometry.A = L*D;
Geometry.V = Geometry.L*Geometry.A;

RconvFloorEnclosure = conductionHorizClosedEnclosures(Tinv,Tout,Geometry,[]);

qconv = (Tinv - Tout) / RconvFloorEnclosure;

disp(['Convective Heat : ' num2str(qconv/1000) ' kW'])

Geometry.L = L;
Geometry.Di = D;
Geometry.b = 0.02;
Geometry.theta = pi; %semicylinder
RconvCylinderEnclosure = naturalConvectionCylinderClosedEnclosures(Tinv, Tout, Geometry, []);


qconv = (Tinv - Tout) / RconvCylinderEnclosure;

disp(['Convective Heat : ' num2str(qconv/1000) ' kW'])


Tfloor = linspace(0, 25, 1);
Tout = 10; % el aire afuera esta a 0 (solo para probar la funcion)
% hay que revisar el modelo de Tout. no sabemos 100% si podemos tomar como
% Tsky la temperatura de afuera


Geometry.A = L*D;

MaterialProperties.emissivity = 0.85; % 

RReradiation = reradiation(Tfloor, Tout, Geometry, MaterialProperties);

qReradiation = (Tfloor - Tout)/RReradiation;

disp(['Leaving Radiative Heat : ' num2str(qReradiation/1000) ' kW'])

Geometry.L = L;
Geometry.D = D;
Geometry.A = Geometry.L*Geometry.D;
Geometry.P = 2*Geometry.L + 2*Geometry.D;

RconvFloorFlatPlate = naturalConvectionHorizFlatPlate(Tinv, Tout, Geometry, []);

qconv = (Tinv - Tout) / RconvFloorFlatPlate;

disp(['Convective heat : ' num2str(qconv/1000) ' kW'])

% perdidas por ventilacion

V = pi*D^2/4/2*L; %m3

airFlow = V/3/60; %m3/s

Geometry.V = V;
Geometry.airFlow = airFlow;

RVentilation = ventilationLoss(Tinv, Tout, Geometry, []);

qVentilation = (Tinv - Tout)/RVentilation;

disp(['Ventilation heat : ' num2str(qVentilation/1000) ' kW'])



% para terminar el balance falta conectar
% mejor entre los nodos
