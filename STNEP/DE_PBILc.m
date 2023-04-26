function [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]= DE_PBILc(...
       MaxIt,Max_FES,SS,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
       incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
       landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion,varargin)

 Acc = 1e-3; % Precision Deseada
 RunStats = [ ]; % Matriz de estadísticas promediadas
 
%% INICIALIZACIÓN DE LOS PARÁMETROS  DEL "PBILc"
 mejores_individuos=round((SS)/2); % Porcentaje del total de mejores individuos
 Tasa_de_aprendi=0.005;              % Tasa de aprendizaje
%% INICIALIZACIÓN DE LOS PARÁMETROS DEL ALGORITMO "DE"
F1=1;    % Factor de escala F
Cr_rand=0.3;

%% INICIALIZACIÓN DE LOS PARÁMETROS (combinar DE y PBILc)
DE_PBILc=0.9;% valores para seleccionar entre De y PBILc
%% INICIALIZACIÓN DE LOS PARÁMETROS DE DOBLE MUTACIÓN(DE)
% valores para la auto mutación (mutación doble del DE)
auto_mutacion_lim_inf=0.3;  %limite superior
auto_mutacion_lim_sup=0.3;  %limite inferior

 %% INICIALIZACIÓN DE VECTORES Y MATRICES NECESARIAS EN EL ALGORITMO        
     xopt=zeros(dim,nexp); % Matriz que almacena la mejor topología en cada iteración
     FOPT=zeros(MaxIt,nexp);      % Matriz que almacena el costo de la mejor topologia en cada iteración
     iteraciones_cada=zeros(MaxIt,nexp); % Matriz para almacenar los el numero de evaluaciones de la funcion objetivo en cada iteracion
     mpc=loadcase(sistema);  % carga algunas variables necesarias del sistema de prueba bajo estudio esta es una función de matpower 
     Comp_Indu_Capa_inicial=mpc.bus(:,6); % lectura para establecer una dimensión de acuerdo a la mpc.bus
     Compensacion_final_suma=ones(size(Comp_Indu_Capa_inicial,1),nexp); % establecer una dimension de acuerdo a la mpc.bus
     fevalcount=zeros(1,nexp);
 %===================================
 % Lazo del número de experimentos
 %===================================
 for T = 1:nexp     
    tic %Inicializar tiempo de simulacion para el ensayo T>1  
    Costos= mpc.ne_branch(:,14); 
 % INICIALIZAR PARAMETROS
    success = 0; % Bandera de éxito
    iter = 1; % Contador de iteraciones
    fevalcount_dina = 0; % Contador del # del Evaluaciones de la Función Objetivo
    STOP = 0; % Bandera de parada del experimento
    contador=0; % Contador de valor mínimo repetido para incorporar el algoritmo del caos
    contador_parada=0; % Inicialización del contador del numero de veces que se repite el valor minimo sin cambiar
     %% Inizialicacion de matrices para el DE           
    U=zeros(dim,SS); %%Incializacion del vector para individuos mutados
    V=zeros(dim,SS); %Inicialización del vector los individuos cruzados         
     
%% Inicialización del vector de probabilidad (PBILc)
    poblacion2=zeros(dim, SS); %Inicialización de la matriz  usada por el PBILc
    Dominio_Maximo=Xmax(:,1);
    Dominio_Minimo=Xmin;
    Rango=Dominio_Maximo-Dominio_Minimo;
    ProbVec =((Dominio_Minimo+2* Rango.*rand(dim, 1))); %vector de probabilidad (media)       
    ProbVec1 =(2*ones(dim,1));      % vector de probabilidades  (varianza)

%% Inicializacion de matrices para almacenar la compensacion shunt
Comp_Indu_Capa_etap=zeros(size(Comp_Indu_Capa_inicial,1),SS);% falta hallar la dimencion de la matriz 
solucion_compe=zeros(size(Comp_Indu_Capa_inicial,1),SS);
f_U=zeros(1,SS);
costo_lineas=zeros(1,SS);
penalizacion=zeros(1,SS);


%% Llamar a función creacioninicial para obtener soluciones iniciales
[poblacion,fevalcount_creacion]=creacioninicial(fhd,...
sistema,Comp_Indu_Capa_etap,Contin_Segun_Dimensi_Final,...
incluir_contingencias,incluir_perdidas,factor_perdidas,...
compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
SS,Xmin,Xmax,dim,opcion); 
    %% Evaluar el flujo optimo de potencia para cada individuos en cada etapa
        Eva=zeros(1,SS); 
        for i=1:SS
            [costo_lineas(i),penalizacion(i),Eva(i),solucion_compe(:,i)]=feval(fhd,i,poblacion(:,i),poblacion,...
            sistema,Comp_Indu_Capa_etap,Contin_Segun_Dimensi_Final,...
            incluir_contingencias,incluir_perdidas,factor_perdidas,...
            compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
            SS);            
        end 

    %% Actualizar la compensación shunt requerida en la etapa Qp para la topología i       
        for i=1:SS
            Comp_Indu_Capa_etap(:,i)=solucion_compe(:,i);
        end
    %% Actualizar el número de evaluaciones de la función objetivo     
    fevalcount_dina=sum(Eva(+1:SS))+fevalcount_dina;    
  
%% Actualizar el valor de costo total (líneas y penalización) para todas las etapas n_anos
    f_costo_X=zeros(1,SS);
    for i=1:SS
        f_costo_X(i)=costo_lineas(i)+penalizacion(i);            
    end 
              
%% ORDENAR LOS INDIVIDUOS EN FUNCIÓN DEL COSTO   
    [fxopt, posicion_mejor] = sort(f_costo_X,'ascend');% Ordenar individuos en base a la funcion objetivo 
    bueno=fxopt(1);% Mejor valor de la funcion objetivo
    FOPT(iter,T)=bueno;% Guardar mejor valor de la funcion objetivo de cada iteracion de cada test
 
%% Actualización de poblacion2 usado para actualizar el vector de probabilidades (PBILc)                
    for i=1:SS
        poblacion2(:,i)=poblacion(:,posicion_mejor(i));     
    end
         
%% Actualización de la mejor topología y su compensación 
    xopt(:,T)=poblacion2(:,1); %topologías que proporciona la mejor solucion para el estado dinamico       
    Compensacion_final_suma(:,T)=Comp_Indu_Capa_etap(:,posicion_mejor(1));  
%% Actualizar el Contador del numero de evaluaciones de la funcion objetivo hasta ahora
   iteraciones_cada(iter,T)=fevalcount_dina;

%======================
% LAZO DE EVOLUCIÓN
%======================
 while (STOP == 0)    
%% ACTUALIZAR CONTADOR DE ITERACIONES

    iter = iter+1;        
%% Actualizar la mejor solución anterior (valor mínimo anterior)         
    cost_ant=bueno;    
%% Actualización del vector de probabilidad (PBILc)
        ProbVec(:,1)=(1-Tasa_de_aprendi)*ProbVec(:,1)+Tasa_de_aprendi*(poblacion2(:,1)+poblacion2(:,2)-poblacion2(:,SS));
        ProbVec1(:,1)=(1-Tasa_de_aprendi)*ProbVec1(:,1)+Tasa_de_aprendi *std(poblacion2(:,1:mejores_individuos),1,2);
 
%% VERIFICAR QUE EL VECTOR DE PROBABILIDAD NO CONVERJA (PBILc)    
    for i=1:dim
        ProbVec(i,1) = max(min(ProbVec(i,1),Dominio_Maximo(i,1)),Dominio_Minimo(i,1));
    end 
                         
%% GENERAR UN NUEVO DESCENDIENTE
      for i=1:(SS)
          if rand()<DE_PBILc % Utilizar mutacion del "DE"
             indices=randperm(SS,4);
             R1=indices(1); % Obtener numero aleatorios
             R2=indices(2); % Obtener numero aleatorios
             R3=indices(3); % Obtener numero aleatorios
             R4=indices(4);
 %% CREACIÓN DEL VECTOR DONADOR
            mutacion_aleatoria=rand(); % Generar un numero aleatorio
            auto_mutacion=(auto_mutacion_lim_inf+(auto_mutacion_lim_sup-auto_mutacion_lim_inf)*(rand));
            if mutacion_aleatoria<auto_mutacion % Utilizar mutacion trigonometrica
                pert_total=f_costo_X(R1)+f_costo_X(R2)+f_costo_X(R3);
                pert1=f_costo_X(R1)/pert_total;
                pert2=f_costo_X(R2)/pert_total;
                pert3=f_costo_X(R3)/pert_total;          
%% Generar un individuo usnado mutacion trigonometrica
                   for k=1:dim 
                        V(k,i)=round((poblacion(k,R1)+poblacion(k,R2)+poblacion(k,R3))/3+...
                                              (pert2-pert1)*(poblacion(k,R1)-poblacion(k,R2))+...
                                              (pert3-pert2)*(poblacion(k,R2)-poblacion(k,R3))+...
                                              (pert1-pert3)*(poblacion(k,R3)-poblacion(k,R1)));  
                   end
%% Verificar los límites de cada individuos en cada etapa 
                    for k=1:dim
                        if V(k,i)<Xmin(k,1)
                            V(k,i)=Xmin(k,1);                      
                        end 
                        if V(k,i)>Xmax(k,i)
                            V(k,i)=Xmax(k,i);                      
                        end                                
                    end                                                         
     
%% UTILIZAR CRUZAMIENTO BINOMIAL
                   U(:,i)=poblacion(:,i);
                   jrand=randi(dim);
                   for k=1:dim
                      if (rand<=Cr_rand)||(k==jrand)
                          U(k,i)=V(k,i);
                      end
                   end              
                
            else                     
%% Generar un individuo usando mutacion difenrecial (DE)
                   for k=1:dim 
                           V(k,i)=round(poblacion(k,posicion_mejor(1))+F1*(poblacion(k,R1)+...
                                           poblacion(k,R2)-poblacion(k,R3)-poblacion(k,R4)));  
                   end 

%% Verificar que los individuos no sobrepasen los limites de su estados        
                    for k=1:dim
                        if V(k,i)<Xmin(k,1)
                           V(k,i)=Xmin(k,1);                      
                        end 
                        if V(k,i)>Xmax(k,i)
                           V(k,i)=Xmax(k,i);                      
                        end                        
                    end                                                         
      
%% UTILIZAR CRUZAMIENTO BINOMIAL
               U(:,i)=poblacion(:,i);
               jrand=randi(dim);
                   for k=1:dim
                      if (rand<=Cr_rand)||(k==jrand)
                          U(k,i)=V(k,i);
                      end
                   end               
                               
            end         
%% VERIFIAR QUE EL VECTOR DONADOR NO SUPERE LOS LIMITES DE BUSQUEDA                                                 
          else   % GENERAR LA NUEVA POBLACION UTILIZANDO PBILc                                
                U(:,i)=round(normrnd(ProbVec(:,1),ProbVec1(:,1)));                       
                               
          end
%% VERIFIAR QUE EL VECTOR DONADOR NO SUPERE LOS LIMITES DE BUSQUEDA                       
                    for k=1:dim
                        if U(k,i)<Xmin(k,1)
                           U(k,i)=Xmin(k,1);                      
                        end 
                        if U(k,i)>Xmax(k,i)
                           U(k,i)=Xmax(k,i);                      
                        end                         
                    end                                                                      
      end 
%% Determinar si el indiivduo actual presenta mejorvalor por adicion de lineas que el individuo anterior      
   ind_eval=zeros(1,SS);
   costo_temp=zeros(1,SS);
          for i=1:SS
              [costo_temp(i)]=costo_adicion_MT(U(:,i),sistema,dim,Costos);
              if costo_temp(i)<f_costo_X(i)
                  ind_eval(i)=i;
              else
                  ind_eval(i)=0;
              end
              f_U(i)=f_costo_X(i);        
          end      
      Comp_Indu_Capa_etap1=Comp_Indu_Capa_etap; 
      Eva=zeros(1,SS);  
              
%% Evaluar el flujo optimo de potencia para cada individuos en cada etapa
    parfor i=1:SS
        if i==ind_eval(i) 
            [costo_lineas(i), penalizacion(i),Eva(i),solucion_compe(:,i)]=feval(fhd,i,U(:,i),U,...
            sistema,Comp_Indu_Capa_etap,Contin_Segun_Dimensi_Final,...
            incluir_contingencias,incluir_perdidas,factor_perdidas,...
            compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,SS) 
        end                                                         
    end    
    %% Actualizar la compensación shunt requerida en la etapa Qp para la topología i       
        for i=1:SS
            if i==ind_eval(i) 
                Comp_Indu_Capa_etap(:,i)=solucion_compe(:,i);
                f_U(i)=costo_lineas(i)+penalizacion(i);            
            end             
        end
%% Actualizar el numero de evaluaciones de la funcion objetivo      
          fevalcount_dina=sum(Eva(1:SS))+fevalcount_dina;             
        
%% Aplicar selecion (DE)         
    for i=1:SS   
        if f_U(i)<f_costo_X(i)
            f_costo_X(i)=f_U(i);
            poblacion(:,i)=U(:,i);
        else
            Comp_Indu_Capa_etap(:,i)=Comp_Indu_Capa_etap1(:,i);            
        end
    end  
 %%  ORDENAR LOS INDIVIDUOS EN FUNCION DEL COSTO       
    [fxopt, posicion_mejor] = sort(f_costo_X,'ascend');
    bueno=fxopt(1);
%% Ordenar a la poblacion actual para actualizar el espacio de creencias (PBILc)
     for i=1:SS
         poblacion2(:,i)=poblacion(:,posicion_mejor(i));     
     end
%% Actualizacion de la mejor topologia y su compensacion 
     xopt(:,(T))=poblacion2(:,1); %topologias que proporciona la mejor solucion para el estado dinamico       
     Compensacion_final_suma(:,(T))=Comp_Indu_Capa_etap(:,posicion_mejor(1));  
                       
%% Guardar la mejor solucion  
   FOPT(iter,T)=bueno;
   figure(1);
   plot(FOPT(:, T)');
   xlabel({['Iter: ',num2str(iter)],['Experimento: ',num2str(T)],['Valor F.O: ',num2str(fxopt(1))]});
   ylabel('Fitness (x 1000 US$)');
   hold on;
%% CONTADOR DEL CRITERIO DE PARADA
    if cost_ant==FOPT(iter, T)
        contador_parada=1+contador_parada;
    else
        contador_parada=1;
    end 
%% Actualizacion de la variable contador 
    if cost_ant==FOPT(iter, T)
        contador=1+contador;
    else
        contador=1;
    end    
     iteraciones_cada(iter,T)=fevalcount_dina;  
%% VERIFICAR CRITRERIO DE PARADA               
     if ((fevalcount_dina >= Max_FES)||(iter >= MaxIt)||(bueno <= Valor_min_referencia+Acc)||contador_parada>Numero_iter_sin_cambiar)
         STOP = 1;
        if (bueno <= Valor_min_referencia+Acc)
            success = 1;
        end
     end
 end % FINALIZAR EVOLUCION 

    fevalcount(1,T)=fevalcount_dina+fevalcount_creacion;
    C_FINAL_CADA_ETAPA=Compensacion_final_suma;
    Tiempo_total_de_simlacion=toc;               %tiempo total de simulacion
    
    %% GUARDAS RESULTADOS EN WIZARD(MATLAB)
    [RunStats]=Guardar_etapa(RunStats, xopt,T,sistema,...
    incluir_contingencias,incluir_perdidas,...
    factor_perdidas,compensacion_anual_optim,landa,C_FINAL_CADA_ETAPA,...
    Comp_Indu_Capa_inicial,success,iter,fevalcount,nexp,...
    Contin_Segun_Dimensi_Final,dim,Particle_Number_original,Metaheuristica_Usada,...
    fxopt,FOPT,opcion,iteraciones_cada,...
    Tiempo_total_de_simlacion); 
 end % FINALIZAR LAZO DEL NUMERO DE EXPERIMENTOS  

%%   IMPRIMIR ESTADISTICAS PROMEDIO
 fprintf('\nESTADISTICAS PROMEDIO:\n');
 SucRunStats = RunStats(RunStats(:,1)==1,:);
 fprintf('EXITOS = %3d/%3d\n', sum(RunStats(:,1)),nexp);
 fprintf('ITERACIONES PROMEDIO = %8.2f\n', mean(SucRunStats(:,2)));
 fprintf('ITERACIONES DESVIACION St= %8.2f\n', std(SucRunStats(:,2)));
 fprintf('EVALUACIONES PROMEDIO DE FUNCION OBJETIVO . = %8.2f\n', mean(SucRunStats(:,3)));
 fprintf('DESVIACION St DE EVALUACION DE FUNCION OBJETIVO = %8.2f\n', std(SucRunStats(:,3)));
 fprintf('VALOR OPTIMO= %8.2f\n', fxopt(1));
 fprintf('EVALUACIONES DE LA FUNCION OBJETIVO DEL ULTIMO EXPERIMENTO= %8.2f\n', fevalcount);

 
end
 