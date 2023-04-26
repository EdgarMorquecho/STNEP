function [resistencia_total, reactancia_total, susceptancia_total, limite_flujo]=Topologia_Red_aux(poblacion, sistema, Qp)

%% Adecuación de los parámetros para el cálculo de las restricciones
    
    x_posicion=poblacion(:,Qp);

    if Qp>1
        LineasIniciales=(poblacion(:,Qp-1))'; % recoje los costos y la topologia inicial de la etapa anterior
    else 
      [LineasIniciales,~]=TestSystems(sistema); % recoje los costos y la topologia inicial 
    end
    
    define_constants; %comodin para acceder a datos
    mpc=loadcase(sistema); %cargar el sistema a resolver
    
    numero_ramas=numel(LineasIniciales);
    LT_nuevas=x_posicion-LineasIniciales';
    Lineas_Totales=LT_nuevas+LineasIniciales'; %lineas totales=x_posicion
  
    resistencia_lineas=mpc.branch(:,3);
    reactancia_lineas=mpc.branch(:,4);
    susceptancia_paralelo=mpc.branch(:,5);
    
    susceptancia_total=zeros(numero_ramas,1);
    reactancia_total=zeros(numero_ramas,1);
    resistencia_total=zeros(numero_ramas,1);

for i=1:numero_ramas
    %%% Adecuacion de susceptancia en paralelo
    susceptancia_total(i,1)=susceptancia_paralelo(i,1)*Lineas_Totales(i,1);
    %%% Adecuacion de resistencia de linea
    if resistencia_lineas(i,1)~=0
        resistencia_total(i,1)=resistencia_lineas(i,1)/Lineas_Totales(i,1);
    else
        resistencia_total(i,1)=0;
    end
    %%% Adecuacion de reactancia de linea
    if reactancia_lineas(i,1)~=0
        reactancia_total(i,1)=reactancia_lineas(i,1)/Lineas_Totales(i,1);
    else
        reactancia_total(i,1)=0;
    end
    
    if Lineas_Totales(i,1)==0
        reactancia_total(i,1)=0;
        resistencia_total(i,1)=0;
        %nuevo sin considerar anteriormente
        %%susceptancia_total(i,1)=0;
        %
    end
end

%%% Adecuacion de limite de flujos
    limite_inicial=mpc.branch(:,6);             %%%% Limite MVA PARA MODELOS AC

%limite_flujo=x_posicion(1:numero_ramas,:).*limite_inicial;   %multiplica el numero de lineas de la topologia actual por su límite de potencia, dando el limite de flujo total por derecho de transmisión
limite_flujo=x_posicion(:,1).*limite_inicial;    %multiplica el numero de lineas de la topologia actual por su límite de potencia, dando el limite de flujo total por derecho de transmisión 