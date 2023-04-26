function [RunStats]=Guardar_etapa(RunStats, xopt,T,sistema,incluir_contingencias,incluir_perdidas,...
factor_perdidas,compensacion_anual_optim,landa,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,success,iter,fevalcount,~,...
Contin_Segun_Dimensi_Final,dim,Particle_Number_original,Metaheuristica_Usada,fxopt,FOPT,opcion,iteraciones_cada,...
Tiempo_total_de_simlacion)  % Graba el archivo .mat con el resultado de algunas variables)
    RunStats(T,:) = [success iter fevalcount(1,T) Tiempo_total_de_simlacion];
    SucRunStats = RunStats(RunStats(:,1)==1,:);
    ITER_PROMEDIO = mean(SucRunStats(:,2));
    ITER_DESV_St=std(SucRunStats(:,2));
    EVA_PROMEDIO_FO= mean(SucRunStats(:,3));
    DESV_St__EVA_FO= std(SucRunStats(:,3)); 
      
%[LineasIniciales,Costos]=TestSystems(sistema); % Algunas variables necesarias que se necesitan almacenar     
Comp_final_guardado_etapa=sum(C_FINAL_CADA_ETAPA);

%% CALCULO DEL NUMERO DE LINEAS ADICIONADAS EN CADA PERIODO DINAMICO   
Lineas_adicionadas_cada_estado1=xopt;
mpc=loadcase(sistema);
   Costo_Lineas_cada_estado1=zeros(dim,T);
    for i=1:T
        Lineas_adicionadas_cada_estado1(:,i)=xopt(:,1);
        Costo_Lineas_cada_estado1(:,i)=mpc.ne_branch(:,14).*Lineas_adicionadas_cada_estado1(:,i); 
    end
%----------------------------------
mpc=loadcase(sistema); 
Barras_branch_col=mpc.ne_branch(:,1:2);
Lineas_adicionadas_cada_estado=[Barras_branch_col Lineas_adicionadas_cada_estado1];
Costo_Lineas_cada_estado=[Barras_branch_col Costo_Lineas_cada_estado1];
Barra_generacion_1=mpc.bus(:,1);
%-----------------------------------------------------
%% Calculo del costo de lineas adicionadas en cada periodo dinamico
Costo_total_lineas_cada_estado=sum(Costo_Lineas_cada_estado1); %Calculo del costo total de lineas para cada periodo
Costo_Total_lineas_todos_estados=zeros(1,T);
    for i=1:T
        Costo_Total_lineas_todos_estados(i)=sum(Costo_total_lineas_cada_estado(T:T)); % calculo del costo por adicion para todos los estados dinamicos
    end
%% CALCULO DEL COSTO DE OPERACION PARA CADA ESTADO DINAMICO
Perdidas_Lineas_cada_estado1=zeros(dim,T);% variable  para guardar los resultados de perdidas en el sistema
status_linea=zeros(dim,T);% variable para ver que lineas fueron adcionadas al final del proceso de optimizacion
Costo_Gen_Pot_Activa_etapa=zeros(1,T);
Costo_Gen_Pot_Activa_todas_etapas=zeros(1,T);
Gen_Pot_Activa_etapa=zeros(1,T);
Gen_Pot_Activa_todas_etapas=zeros(1,T);
     mpc=loadcase(sistema);  % carga algunas variables necesarias del sistema de prueba bajo estudio esta es una función de matpower              
     Compensacion_actual_verificada=zeros(size(mpc.bus,1),T); %variable para guardar la compensacion para el siguietne a?o
     Comp_shunt_en_cada_barra1=zeros(size(mpc.bus,1),T);
     Costo_individual_compe=zeros(size(mpc.gen,1),T);
     Barra_con_shunt=zeros(size(mpc.gen,1),T);
for N_exp=1:T
        [Gen_Pot_Activa_etapa(N_exp),Costo_Gen_Pot_Activa_etapa(N_exp),...
        Costo_individual_compe(:,N_exp),...
        Barra_con_shunt(:,N_exp),Comp_shunt_en_cada_barra1(:,N_exp)]=Calculo_Comp_Gen_Perd(xopt(:,N_exp),...
        sistema,xopt,compensacion_anual_optim,Compensacion_actual_verificada,...
        Comp_Indu_Capa_inicial,N_exp,incluir_contingencias,Contin_Segun_Dimensi_Final);  
        Compensacion_actual_verificada(:,N_exp)=Comp_shunt_en_cada_barra1(:,N_exp);     
    Gen_Pot_Activa_todas_etapas(N_exp)=sum(Gen_Pot_Activa_etapa(N_exp:N_exp)); % Calcula la generacion total de todos los estados    
    Costo_Gen_Pot_Activa_todas_etapas(N_exp)=sum(Costo_Gen_Pot_Activa_etapa(N_exp:N_exp));
end
       
%% CALCULO DE LAS PERDIDAS (para cada estado)   
         for i=1:dim
             for j=1:T
                 if  status_linea(i,j)==0
                     Perdidas_Lineas_cada_estado1(i,j)=0;
                 end
             end
         end          
%% Calculo de las perdidas totales y el costo para cada numero de experimento (nexp)
Perdidas_total_cada_estado=sum(Perdidas_Lineas_cada_estado1); % Sumar las perdidas de cada derecho de trasmision (Perdidas totales)
Perdidas_Lineas_cada_estado=[Barras_branch_col Perdidas_Lineas_cada_estado1];
Perdidas_total_todos_estados=zeros(1,T);
Cost_Perdidas_total_todos_estados=zeros(1,T);
Costo_Perdida_cada_estado=zeros(1,T);
for i=1:T
    Perdidas_total_todos_estados(i)=sum(Perdidas_total_cada_estado(i:i));
      if incluir_perdidas==1
             Costo_Perdida_cada_estado(i)=Perdidas_total_cada_estado(i)*landa*factor_perdidas*8760;

      else
             Costo_Perdida_cada_estado(i)=0;
            
      end
      Cost_Perdidas_total_todos_estados(i)=sum(Costo_Perdida_cada_estado(i:i));
end

%% CALCULO DE LA COMPENSACION DE REACTIVOS PARA UNA TOPOLOGIA ESPECIFICA EN CADA ESTADO 
     mpc=loadcase(sistema);  % carga algunas variables necesarias del sistema de prueba bajo estudio esta es una función de matpower              
     Comp_Indu_Capa_barra=zeros(size(mpc.bus,1),T); %variable para guardar la compensacion para el siguietne a?o
     Comp_Inductor=zeros(size(mpc.bus,1),T); %variable para guardar la compensacion para el siguietne a?o
     Comp_Capacitor=zeros(size(mpc.bus,1),T); %variable para guardar la compensacion para el siguietne a?o
     Costo_shunt_total_barra_result=zeros(size(mpc.bus,1),T); %variable para guardar el costo de la compensacion para el siguietne a?o
    [Posicion_Nodos_unicos,~]=unique(mpc.gen(:,1),'rows');%Determina barras no repetidas
     for i=1:T
            for j=1:size(mpc.bus)
                Comp_ficticia=find(mpc.gen(:,1)==Posicion_Nodos_unicos(j));
                if size(Comp_ficticia,1)==1
                    Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i)=0; 
                elseif size(Comp_ficticia,1)==2
                    Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i)=(Costo_individual_compe(Comp_ficticia(1),i)+Costo_individual_compe(Comp_ficticia(2),i));                    
                    if Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i)>=0
                      Comp_Capacitor(Posicion_Nodos_unicos(j),i)= Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i);
                      Costo_shunt_total_barra_result(Posicion_Nodos_unicos(j),i)=Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i)*Costo_individual_compe(Comp_ficticia(1),1);  
                    else 
                      Comp_Inductor(Posicion_Nodos_unicos(j),i)= Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i); 
                      Costo_shunt_total_barra_result(Posicion_Nodos_unicos(j),i)=Comp_Indu_Capa_barra(Posicion_Nodos_unicos(j),i)*Costo_individual_compe(Comp_ficticia(2),1);    
                    end                          
                end
            end           
     end 

  Comp_Inductor_Total=sum(Comp_Inductor);% Compensacion shunt inductiva para cada estado
  Comp_Capacitor_Total=sum(Comp_Capacitor);% Compensacion shunt capacitiva para cada estado
  Comp_shunt_costo_etapa=sum(Costo_shunt_total_barra_result);% Compensacion shunt total para cada estado
  Comp_shunt_costo_total=zeros(1,T); 
  Compensacion_shunt_total_Stage=zeros(1,T); 
  Compensacion_shunt_total=sum(Comp_shunt_en_cada_barra1);% Compensacion shunt total para cada estado  
  for i=1:T
      Comp_shunt_costo_total(1,i)=sum(Comp_shunt_costo_etapa(1,i:i));
      Compensacion_shunt_total_Stage(i)=sum(Compensacion_shunt_total(i:i));      
  end
  Comp_shunt_en_cada_barra=[Barra_generacion_1 Comp_shunt_en_cada_barra1];   





    %% VARIABLES NECESARIAS PARA GUARDAR EN UN DOCUMENTO .MAT  
    Particle_name=num2str(Particle_Number_original);  % Transforma el nUmero de individuos en caracteres
    PERDAS=num2str(incluir_perdidas);  % Transforma el n?mero de individuos en caracteres
    SI_continge=num2str(incluir_contingencias);  % Transforma el n?mero de individuos en caracteres
    E_dinamicas=num2str(1);  
    Test_numero=num2str(T);     
%% HORA EN LA QUE FINALIZA LA SIMULACION     
    Hora_Actual=datestr(now);
    Hora_Actual1=Hora_Actual(13:14);
    Hora_Actual2=Hora_Actual(16:17);
    Hora_Actual3=Hora_Actual(19:20);

%% NOMBRE DEL ARCHIVO A GUARDAR    
filename=strcat(sistema,'_',Particle_name,'_',date,'_',Metaheuristica_Usada,'_',SI_continge,'_',PERDAS','_','Hora','_',Hora_Actual1,'_',Hora_Actual2,'_',Hora_Actual3,'_','N_Etap','_',E_dinamicas,'_',Test_numero);  % Escribe un archivo .mat    
save(filename,'Comp_Capacitor',...
                'Comp_Inductor','iteraciones_cada','opcion','Lineas_adicionadas_cada_estado',...
                'Costo_Lineas_cada_estado','Costo_total_lineas_cada_estado',...
                'Costo_Total_lineas_todos_estados','Perdidas_Lineas_cada_estado','Perdidas_total_cada_estado',...
                'Perdidas_total_todos_estados','Costo_Perdida_cada_estado','Cost_Perdidas_total_todos_estados',...
                'Compensacion_shunt_total_Stage','Comp_Inductor_Total','Comp_Capacitor_Total','Comp_final_guardado_etapa',...
                'Compensacion_shunt_total','Comp_shunt_costo_total','Comp_shunt_en_cada_barra','Gen_Pot_Activa_todas_etapas',...
                'Costo_Gen_Pot_Activa_todas_etapas','Costo_Gen_Pot_Activa_etapa','Gen_Pot_Activa_etapa','xopt',...
                'C_FINAL_CADA_ETAPA','landa', 'fxopt',...
                'fevalcount', 'FOPT','RunStats','ITER_PROMEDIO','ITER_DESV_St','EVA_PROMEDIO_FO',...
                'DESV_St__EVA_FO') % Graba el archivo .mat con el resultado de algunas variables) % Graba el archivo .mat con el resultado de algunas variables
            
% Lineas_adicionadas_cada_estado: guarda las lineas adicionadas en cada estado
% Costo_Lineas_cada_estado: Calcula el costo de lineas por cada linea adicionada en cada derecho de trasmision 
% Costo_total_lineas_cada_estado: Calcula el costo total de lineas adicionadas por cada estado
% Costo_Total_lineas_todos_estados: Calcula el costo total por lineas adicionadas durante todas las etaps dinamicas
% Costo_Lineas_cada_estado: calcula el costo poradicion de lineas para cada derecho de trasnmision    
% Costo_lineas_cada_estado: calcula el costo total por adicion de lineas para cada estado  
% Costo_Total_lineas_todos_estados: calcula el costo total por adicion de lineas para todas las etapas 
% Perdidas_Lineas_cada_estado: Calcula las perdidas en cada derecho de trasnmision para cada etapa
% Perdidas_total_cada_estado: Calcula las perdidas totales para cada estado 
% Perdidas_total_todos_estados: Calcula las perdidas totales para todos los estados 
% Costo_Perdida_cada_estado: Calcula el csoto por perdidas para cada estado o etapa
% Cost_Perdidas_total_todos_estados: Calcula las perdidas totales para cada ensayo (nexp)  
% Comp_Inductor_Total: Calcula la compensacion shunt inductiva para cada estado
% Comp_Capacitor_Total: Calcula la compensacion shunt capacitiva para cada estado    
% Comp_shunt_costo_etapa: Calcula el costo por compensacion shunt para cada etapa o estado     
% Comp_shunt_costo_total: Calcula el costo total de compensacion para cada ensayo (nexp)     
% Compensacion_shunt_total_Stage: Calcula la compensacion shunt total para cada ensayo (nexp)
% Compensacion_shunt_total: Calcula la compensacion shunt para cada etapa o estado
% Comp_shunt_en_cada_barra: Compensacion shunt tanto inductiva compo capacitica para cada estado  
% Gen_Pot_Activa_todas_etapa:Calcula la generacion activa total para todas las n etapas
% Costo_Gen_Pot_Activa_todas_etapas:Calcula el costo por generacion activa total para todas las n etapas 
% Costo_Gen_Pot_Activa_etapa: Calcula el costo de generacion activa total para todas cada etapas
% Gen_Pot_Activa_etapa:Calcula la generacion activa total para cada etapas
% Comp_final_guardado_etapa: calcula la compensacion shunt total de cada etapa(este debe ser similar a Comp_shunt_costo_etapa)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
