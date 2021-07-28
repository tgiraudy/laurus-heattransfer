%% Ubicacion

Options.Ubicacion = [-34.4 -58.6];
Options.GMT = -3;
Options.Albedo = 1;

%% elements

elements = [1 2   % Conveccion natural interna
            2 3   % Conveccion en enclosure
            1 4]; % Conveccion natural con piso


        
nEle = size(elements, 1);
nNodEle = size(elements, 2);
nNod = max(max(elements));
nDofNod = 1;
nDofEle = nDofNod*nNodEle;
nDofTot = nDofNod*nNod;
dofs = reshape(1:nDofTot, nDofNod, [])';


Kele = zeros(nDofEle, nDofEle, nEle);
T = [25 22 5.5 15]';
prevT = T;
q = zeros(nDofTot, 1);
interstepTol = 1e-3;
interstepError = 1;

beta = 1;
fixed = false(nDofTot, 1);
fixed([4 6 7]) = true;

startDay = 180;
nDays = 3;



endTime = 24*nDays; %h
time = 0:endTime;
nSteps = length(time);

TinfData = 5*ones(nSteps, 1); % provisional
T0Data = 10*ones(nSteps, 1); 
radiationVector = zeros(1, nSteps);
Tinv = zeros(nSteps, 1);
Tinv(1) = 20;



%% branches al invernadero

[branchElements, branchDirection] = find(elements == 1); % 1 es el indice del aire ambiente del invernadero


%% Geometria        
L = 6;
D = 4;
V = L*pi*D^2/4/2;

bAislacion = 0.02;
bSuelo = 0.2;
        
c2k = @(T) T + 273.15;        
% si se agregan costados tiene que entrar en el calculo de elementos 1, 2 y BC;

%% Propiedades

solarAbsorptivity = 0.5;
floorEmissivity = 0.5;
air = Air();

airMass = air.rho*V;


T(1) = Tinv(1);

T = c2k(T);

for tIndex = 1:nSteps; %horas
    
    day = startDay + floor(time(tIndex)/24);
    hour = rem(time(tIndex),24);
    %
    timeStep = 1;
    radiationVector(tIndex) = SunRadiation([0 1 0], hour, day, Options);
    Tinf = TinfData(tIndex);
    T0 = T0Data(tIndex);
    %         hr = hrData(tIndex);
    Tdp = Tinf; % Temperatura dew point (EN REALIDAD TIENE QUE SER FUNCION DE Tinf y hr);
    epsilonSky = 0.711 + 0.56*(Tdp/100) +  0.73*(Tdp^2/100)^2; % berdahl and martin 1984 REVISAR
    Tsky = (epsilonSky)^.25*Tinf;
    deltaT = zeros(nDofTot, 1);
    
    sunRadiationHeat =  radiationVector(tIndex)*L*D* solarAbsorptivity;
<<<<<<< Updated upstream
    boundaryConditions = [Tinf T0 Tsky sunRadiationHeat];

    [] = newtonRaphson(T,  bc
    
    %% branches
    
    qBranches = 0;
    

    
=======
    boundaryConditions = [c2k([Tinf T0 Tsky]) sunRadiationHeat,airMass];
        
    T = newtonRaphson(T,prevT,boundaryConditions,elements,timeStep);
>>>>>>> Stashed changes
    
    prevT = T; 
    disp(T(1))
end