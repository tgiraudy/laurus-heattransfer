clear
close all
clc

data = [0 ,	28 
        1 ,	28
        2 ,	26 
        3 ,	24 
        4 ,	23 
        5 ,	21 
        6 ,	20 
        7 ,	20 
        8 ,	19 
        9 ,	20 
        10 ,	21 
        11 ,	25 
        12 ,	25 
        13 ,	25 
        14 ,	27 
        15 ,	29 
        16 ,	29 
        17 ,	31 
        18 ,	29 
        19 ,	28 
        20 ,	28 
        21 ,	28 
        22 ,	28 
        23 ,	27];

f2c = @(f)(f-32)*5/9;

data(:,2) = f2c(data(:,2));

time = data(:,1);
T2 = data(:,2); %C

% temperatura menor
T1 = 25; %C

% dimensiones
Di = 10; %m
b = 0.5; %m
L = 15; %m
Do = Di + 2*b;

% propiedades aire
k = 23e-3; % [W/mk]
Pr = 0.71;
gbnu = 1.85e8;

% adimensionales
Gr = gbnu.*(b.^3).*(T1 - T2);
Ra = Gr.*Pr;

%p 321 kreith heat transfer
check = (log(Do./Di)./(b.^.75.*(Do.^(-3/5) + Di.^(-3/5)).^(5/4))).^4.*Ra;
kEff = k*0.386*(log(Do./Di)./(b.^.75.*(Do.^(-3/5) + Di.^(-3/5)).^(5/4))).*(Ra.^(1/4))*((Pr/(0.861+Pr)).^(1/4));

% resistencia medio cilindro
R = log(Do./Di)./(0.5*2*pi*kEff*L);

% calor
q = (T1 - T2)./R;

Energy = trapz(time, q); %kWh
powerAverage = Energy/24; %kW

figure('Color', 'w')
plot(time, q)
grid on
grid minor


floorSurface = L*Di;
floorPerimeter = 2*L+2*Di;
floorEmissivity = 0.84;

TFloor = linspace(0, 25, 20);

characteristicLength = floorSurface/floorPerimeter;

GrL = gbnu.*(characteristicLength.^3).*(T1 - TFloor);
RaL = GrL.*Pr;

NuL = 0.27.*RaL.^(1/4);

hFloor = NuL*k/characteristicLength;
qConvFloor = hFloor.*floorSurface.*(T1 - TFloor);



TOut = 0; % el aire afuera esta a 0 (solo para probar la funcion)

Geometry.A = floorSurface;
MaterialProperties.emissivity = floorEmissivity;

RReradiation = reradiation(TFloor, TOut, Geometry, MaterialProperties);





