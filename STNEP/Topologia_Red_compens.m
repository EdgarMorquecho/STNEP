function [resistencia_total, reactancia_total, susceptancia_total, limite_flujo]=Topologia_Red_compens(x_posicion,sistema,...
    xopt,N_exp)
%% Adecuación de los parámetros para el cálculo de las restricciones



%% Adecuación de los parámetros para el cálculo de las restricciones
numero_ramas=numel(x_posicion);% Determinar la dimension del sistema
mpc=loadcase(sistema); %cargar el sistema a resolver
define_constants; %comodin para acceder a datos
Lineas_Totales=x_posicion;% topologia del sistema
reactancia_lineas=mpc.branch(:,4);
resistencia_lineas=mpc.branch(:,3);
susceptancia_paralelo=mpc.branch(:,5);
susceptancia_total=zeros(numero_ramas,1);
reactancia_total=zeros(numero_ramas,1);
resistencia_total=zeros(numero_ramas,1);
for j=1:numero_ramas
    %% Adecuacion de susceptancia en paralelo
    susceptancia_total(j,1)=susceptancia_paralelo(j,1)*Lineas_Totales(j,1);
    
    %% Adecuacion de reactancia de linea
    if reactancia_lineas(j,1)~=0
    reactancia_total(j,1)=reactancia_lineas(j,1)/Lineas_Totales(j,1);
    else
        reactancia_total(j,1)=0;
    end
    %% Adecuacion de resistencia de linea
    if resistencia_lineas(j,1)~=0
    resistencia_total(j,1)=resistencia_lineas(j,1)/Lineas_Totales(j,1);
    else
        resistencia_total(j,1)=0;
    end
    if Lineas_Totales(j,1)==0
        reactancia_total(j,1)=0;
        resistencia_total(j,1)=0;
    end
end
%% Calculo del limite potencia por derecho de transmision
limite_inicial=mpc.branch(:,6); %
limite_flujo=x_posicion(1:numero_ramas,:).*limite_inicial; %multiplica el numero de lineas de la topologia actual por su límite de potencia, dando el limite de flujo total por derecho de transmisión
