function [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]= TLBO(...
            MaxIt,Max_FES,SS,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
            incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
            landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion,varargin)

 Acc = 1e-3; % Precision Deseada
 RunStats = [ ]; % Matriz de estadísticas promediadas

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
    x_prof=zeros(dim); %Vector que contiene el mejor individuo (profesor)
    x_media=zeros(dim); %Vector que contiene la media de cada etapa 
    poblacion_nueva=zeros(dim,SS);     
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
Eva=zeros(1,SS);
    %% Evaluar el flujo optimo de potencia para cada individuos en cada etapa
        parfor i=1:SS
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
%% Actualización de la mejor topología y su compensación 
    for i=1:SS      
        poblacion_nueva(:,i)=poblacion(:,posicion_mejor(i)); %topologías que proporciona la mejor solucion para el estado dinamico       
    end
%% Actualización de la mejor topología y su compensación 
    xopt(:,T)=poblacion(:,posicion_mejor(1)); %topologías que proporciona la mejor solucion para el estado dinamico       
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
%% Identificar al individuo profesor
        x_prof(:,1) = poblacion(:,posicion_mejor(1)); %Individuo con el mejor valor de la funcion objetivo
        x_media(:,1)=mean(poblacion(:,1:SS),2);
    
%% Fase del profesor
poblacion_nueva=poblacion;
    for i=1:SS
        if i~=posicion_mejor(1)
            TF=round(1+rand);                
            r = rand(dim,1);                
            poblacion_nueva(:,i)=round(poblacion(:,i)+r.*(x_prof(:,1)-TF.*x_media(:,1)));

%% VERIFIAR QUE EL VECTOR DONADOR NO SUPERE LOS LIMITES DE BUSQUEDA                       
            for k=1:dim
                if poblacion_nueva(k,i)<Xmin(k,1)
                   poblacion_nueva(k,i)=Xmin(k,1);                      
                end 
                if poblacion_nueva(k,i)>Xmax(k,i)
                   poblacion_nueva(k,i)=Xmax(k,i);                      
                end                         
            end                                                                                     
        end
    end
%% Determinar si el individuo actual presenta que mejor valor por adicion de lineas que el individuo anterior      
   ind_eval=zeros(1,SS);
   costo_temp=zeros(1,SS);
          for i=1:SS
              [costo_temp(i)]=costo_adicion_MT(poblacion_nueva(:,i),sistema,dim,Costos);
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
            [costo_lineas(i), penalizacion(i),Eva(i),solucion_compe(:,i)]=feval(fhd,i, poblacion_nueva(:,i),poblacion_nueva,...
            sistema,Comp_Indu_Capa_etap,Contin_Segun_Dimensi_Final,...
            incluir_contingencias,incluir_perdidas,factor_perdidas,...
            compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
           SS);  
        end                                                         
    end       
%% Actualizar la compensación shunt requerida en la etapa Qp para la topología i       
        for i=1:SS
            if i==ind_eval(i) 
                Comp_Indu_Capa_etap(:,i)=solucion_compe(:,i);
                f_U(i)=costo_lineas(i)+penalizacion(i);            
            end             
        end  
%% Actualizar el contador del número de evaluaciones de la función objetivo             
              fevalcount_dina=sum(Eva(1:SS))+fevalcount_dina;                    
      
%% Aplicar selección (DE)            
        for i=1:SS    
            if f_U(i)<f_costo_X(i)
                f_costo_X(i)=f_U(i);
                poblacion(:,i)=poblacion_nueva(:,i);           
            else
                Comp_Indu_Capa_etap(:,i)=Comp_Indu_Capa_etap1(:,i);        
            end
        end   
%% Fase del alumno
    for i=1:SS
        ks=randi(SS,1);
        while ks==i
            ks=randi(SS,1);
        end
        if f_costo_X(i)< f_costo_X(ks)
            ro=rand(dim,1);
            poblacion_nueva(:,i)=round(poblacion(:,i)+ro.*(poblacion(:,i)-poblacion(:,ks)));            
        else
            ro=rand(dim,1);
            poblacion_nueva(:,i)=round(poblacion(:,i)+ro.*(poblacion(:,ks)-poblacion(:,i)));            
        end
%% VERIFIAR QUE EL VECTOR DONADOR NO SUPERE LOS LIMITES DE BUSQUEDA                                
        for k=1:dim
            if poblacion_nueva(k,i)<Xmin(k,1)
               poblacion_nueva(k,i)=Xmin(k,1);                      
            end 
            if poblacion_nueva(k,i)>Xmax(k,i)
               poblacion_nueva(k,i)=Xmax(k,i);                      
            end                         
        end       
    end   
%% Determinar si el indiivduo actual presenta mejorvalor por adicion de lineas que el individuo anterior      
   ind_eval=zeros(1,SS);
   costo_temp=zeros(1,SS);
          for i=1:SS
              [costo_temp(i)]=costo_adicion_MT(poblacion_nueva(:,i),sistema,dim,Costos);
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
                [costo_lineas(i), penalizacion(i),Eva(i),solucion_compe(:,i)]=feval(fhd,i,poblacion_nueva(:,i),poblacion_nueva,...
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
            poblacion(:,i)=poblacion_nueva(:,i);
        else
            Comp_Indu_Capa_etap(:,i)=Comp_Indu_Capa_etap1(:,i);            
        end
    end    
                 
 %%  ORDENAR LOS INDIVIDUOS EN FUNCION DEL COSTO       
    [fxopt, posicion_mejor] = sort(f_costo_X,'ascend');
    bueno=fxopt(1);

%% Actualizacion de la mejor topologia y su compensacion 
         xopt(:,T)=poblacion(:,posicion_mejor(1)); %topologias que proporciona la mejor solucion para el estado dinamico       
         Compensacion_final_suma(:,T)=Comp_Indu_Capa_etap(:,posicion_mejor(1));            
              
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
 toc
 
end
 