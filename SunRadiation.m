function [ITot, SunVersor] = SunRadiation(EleNormal,LocalTime,Day,Options)
%Calculo de radiaci�n solar en funcion de la orientacion de una superficie
%plana.
%INPUT: #EleNormal: Array de las normales de los elementos (sistema de coordeanadas x=Este | y=Sky | z=Sur).
%       #LocalTime: Hora(s) del d�a solicitadas.
%       #Day: D�a del a�o.
%       #Options: Struct con los siguientes campos (opcionales).
%           #Ubicaci�n: [Latitud Longitud] (grados) N+ E+.
%           #Albedo: Albedo terrestre.
%           #GMT: Correcci�n al horario local respecto de GMT.
%OUTPUT:#ITot: Irradiaci�n total [xEle, xHora] (W/m2)
%       #SunVersor: Versor radiaci�n solar. (hacia la tierra)
%Referencias: New model to estimate and evaluate the solar radiation - Mghouchi

%Parametros default
% Ubicacacion = [-34.6 -58.3667];         %Coordenadas [Latitud (grados) N+ | Longitud (grados) E+]: BUENOS AIRES                     
Ubicacion = [-33.88 -60.5667];         %Coordenadas [Latitud (grados) N+ | Longitud (grados) E+]: PERGAMINO                    
DTl = -3;                               %Corrimiento hora local GMT: BUENOS AIRES
alb = 0;                                                        %Albedo terrestre

%Chequeo opciones adicionales
if nargin==4 && isstruct(Options)
  if isfield(Options,'Ubicacion'),      Ubicacion = Options.Ubicacion;           end
  if isfield(Options,'Albedo'), alb = Options.Albedo; end
  if isfield(Options,'GMT'), DTl = Options.GMT; end
end

j = Day;                      %D�a del a�o [dias]
Tl = LocalTime;               %Hora del d�a en horario local [horas]

nTl = length(Tl);             %N�mero de horas solicitadas
nEles = size(EleNormal,1);    %N�mero de elementos solicitados

phi = Ubicacion(1);         %Latitud [grados] 
lgt = Ubicacion(2);         %Longitud [grados]


delta = 23.45*sind(0.986*(j+284));   %Declinaci�n solar [grados]
E = 450.8*sin(2*pi*j/365-0.026903) + 595.4*sin(4*pi*j/365+0.352835); %Ec. de tiempo [s]
 
%x=Este | y=Sky | z=Sur
EleIota = acosd(EleNormal(:,2)); %�ngulo respecto de la horizontal
EleGama = atan2d(EleNormal(:,1),-EleNormal(:,3)); %�ngulo gama c/ respecto al norte


iota = repmat(EleIota,1,nTl); %Plano de incidencia: �ngulo respecto de la horizontal 
gama =  repmat(EleGama,1,nTl); %Plano de incidencia: �ngulo respecto del Norte: 0= Norte >0=Hacie el E <0=Hacia O


Tsv = Tl - DTl + (lgt*4 + E/60)/60;                                  %Hora Solar Real [horas]

omega = 15*(12-Tsv);                                             %Angulo horario solar [grados]
 
zen = acosd(sind(phi).*sind(delta)+cosd(phi).*cosd(delta).*cosd(omega));  %Zenith angle
psi = acosd((sind(delta)-cosd(zen).*sind(phi))./(sind(zen).*cosd(phi)));  %Azimut solar 

psi(omega<0) = -psi(omega<0); %Azimut solar: 0= Norte >0=Hacie el E <0=Hacia O

SunVersor = -[sind(zen)'.*sind(psi)', cosd(zen)', -sind(zen)'.*cosd(psi)'];

Io = 1367;                                                      %Constante solar [w/m2]
Ct = 1 + 0.034*cosd(j-2);                                        %Correci�n distnacia tierra-sol []
Tu = 0.796 -0.01 * sind(0.986*(j+284));                           %Tubidez atmosferica []
h = asind(sind(phi).*sind(delta)+cosd(phi).*cosd(delta).*cosd(omega));    %Altura del sol [grados]


Idirecta = Io.*Ct.*Tu.*exp(-0.13./sind(h)).*sind(h); %Radiaci�n directa sup.H [W/m2]
Idifusa = 120.*Tu.*exp(-1./(0.4511+sind(h)));      %Radiaci�n difusa sup.H [W/m2]
Itotal = Idirecta + Idifusa;                    %Radiaci�n total sup.H [w/m2]

Idirecta = repmat(Idirecta,nEles,1);
Idifusa = repmat(Idifusa,nEles,1);
Itotal = repmat(Itotal,nEles,1);
psi = repmat(psi,nEles,1);
h  = repmat(h,nEles,1);

Idirecta(Idirecta<0) = 0;
IdirInc = Idirecta.*((sind(iota).*cosd(psi-gama))./tand(h) + cosd(iota));  %Radiaci�n directa sup.inc [W/m2]

IdirInc(IdirInc<0) = 0;
IdifInc = (1+cosd(iota))./2.*Idifusa+(1-cosd(iota))./2.*alb.*Itotal;         %Radiaci�n difusa sup.inc [W/m2]
ITot = IdirInc+IdifInc; %Radiaci�n total superfice inclinada.

ITot(:,SunVersor(:,2)>0) = 0; %Limpia radiacion de noche
end