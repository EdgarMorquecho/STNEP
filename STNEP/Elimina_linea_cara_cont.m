function [poblacion_anterior1,f_costo_X_anterior1,Comp_Indu_Capa_etap_anterior1,Eva_BL_anterior]...
                =Elimina_linea_cara_cont(poblacion_anterior1,Xmin,Xmax,sistema,...
                incluir_contingencias,incluir_perdidas,factor_perdidas,...
                compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
                SS,dim,...
                Comp_Indu_Capa_etap_anterior1,f_costo_X_anterior1,...
                solucion_compe_anterior1,...
                ~,fhd,Contin_Segun_Dimensi_Final,MapEdgesLocal_2,MapEdgesLocal,Costos)
Eva_BL_anterior=0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ELIMINAR LINEAS BASADO EN EL COSTO (Eliminar una lineas y adicionar una lineas donde el costo de la linea adicionada es menor al costo de la linea eliminada     
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%  Lectura del costo por cada derecho de via de cada etapa              
    C_unitario_linea=Costos;
%% Lectura del mejor individuo actual    
    Cadd=poblacion_anterior1; %Identificar las lineas adicionadas en el individuo para (n_anos=1)                                                                                          
%% Identificar el costo de cada linea adicionada en cada etapa para el individuo actual 
    Costo_linea_reducir=Cadd.*(Costos);
%% Identificar el costo de cada linea en cada etapa para el individuo actual         
        for j=1:dim
            if Costo_linea_reducir(j,1)>0
               Costo_linea_reducir(j,1)= Costos(j,1);
               C_unitario_linea(j,1)=0;%Colocar el costo de cero a las posiciones donde ya se tiene lineas adicionadas
            end
        end        
    indx = find(Costo_linea_reducir~=0); % Identificar posicion donde se adicionan lineas de transmision para el individuo 
    indxr1 = randperm(size(indx,1),round(size(indx,1)*20/100));
    %% Inicio del proceso iterativo                    
    if ~isempty(indxr1)  %Elimina lineas si existen lineas adicionadas en el individuo actual    
%% Iniciar proceso de eliminacion de lineas        
        for indx_j=1:size(indxr1,2)              
%%  Identificar posiciones donde el costo de la linea adicionada  presenta mayor valor que las posciones de las lineas no adicionadas   
            C_unitario_linea_reducido=Costos;%Matriz para identificar posiciones donde el costo de la linea a eliminar no sea mayor al costo de la linea a adicionar
            Costo_maximo_comb=Costo_linea_reducir(indx(indxr1(1,indx_j)), 1);% Costo de la linea a eliminar             
            %C_unitario_linea_reducido(indx(indxr1(1,indx_j)),1)=0;
            for j=1:dim
                if C_unitario_linea_reducido(j,1)>Costo_maximo_comb
                    C_unitario_linea_reducido(j,1)=0;%Colocar el costo de cero a las posiciones donde presenta costo mayor a la combinacion de costos de las lineas adicionadas
                end
            end
%% Identificar posiciones donde aun no se adiciona lineas         
            indx_max = find(C_unitario_linea_reducido~=0); % Identificar posicion donde no se adicionan lineas de transmision para el individuo 
            for indx_k=1:size(indx_max,1)
                particle_mod_anterior = poblacion_anterior1;
%% Eliminar la linea indx_j de lineas adicionadas  
                 particle_mod_anterior(indx(indxr1(1,indx_j)),1)=particle_mod_anterior(indx(indxr1(1,indx_j)),1)-1;%Eliminar una linea 
%% Verificar los limites de la topologia (Limite minimo y maximo)               
                if particle_mod_anterior(indx(indxr1(1,indx_j)),1)<Xmin(indx(indxr1(1,indx_j)),1)
                    particle_mod_anterior(indx(indxr1(1,indx_j)),1)=Xmin(indx(indxr1(1,indx_j)),1);                      
                end                                                                                          
%% Adicionar la  linea indx_k de las posibles lineas a adicionar   
                 particle_mod_anterior(indx_max(indx_k),1)=particle_mod_anterior(indx_max(indx_k),1)+1;%Adicionar una linea 
                      
%% Verificar los limites de la topologia (Limite minimo y maximo)
                if particle_mod_anterior(indx_max(indx_k),1)>Xmax(indx_max(indx_k),1)
                    particle_mod_anterior(indx_max(indx_k),1)=Xmax(indx_max(indx_k),1);                      
                end
                Edges_anterior = cell(1,4);
                Eva=0;  
                Comp_Indu_Capa_etap1_anterior=Comp_Indu_Capa_etap_anterior1;                
                        %% Lineas de trasnmision adicionadas
                        LT_nuevas=particle_mod_anterior; % recoje los costos y la 
                        %% Calcular el costo por lineas adicionadas
                        costo_temp_temp=sum(Costos.*LT_nuevas);
                            if costo_temp_temp<f_costo_X_anterior1
                               ind_eval_temp=1;
                            else
                                ind_eval_temp=0;
                            end
                       
    %% Identificar individuos repetidos y evaluar individuos no repetidos
                if ind_eval_temp==1 
                    try                 
                     key=sprintf('%d;', particle_mod_anterior);
                     sprintf('key %s\n', key);
                     value = MapEdgesLocal(key);
                     value2 = MapEdgesLocal_2(key);            
                     costo_lineas_anterior=value(1);
                     penalizacion_anterior=value(2);
                     solucion_compe_anterior1=value2(:,1);
                     Eva=0; 
                    catch ME
                    if (strcmp(ME.identifier,'MATLAB:Containers:Map:NoKey') || strcmp(ME.identifier,'MATLAB:Containers:TypeMismatch'))      
                       Edges_anterior(1,1)={key}; 
                    else
                        fprintf('Que paso? %s\n', ME.identifier);
                    end
                    end
                end

    %% Evaluar el flujo optimo de potencia para cada individuos en cada etapa
                if ind_eval_temp==1 
                    if (~isempty(Edges_anterior{1,1}))
                        [costo_lineas_anterior, penalizacion_anterior,Eva,solucion_compe_anterior1]=feval(fhd,1,particle_mod_anterior(:,1),particle_mod_anterior,...
                        sistema,Comp_Indu_Capa_etap_anterior1,Contin_Segun_Dimensi_Final,...
                        incluir_contingencias,incluir_perdidas,factor_perdidas,...
                        compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,SS,Xmin);      
                    end 
                end
    %% Actualizar la compensación shunt requerida en la etapa Qp para la topología i  
                if ind_eval_temp==1     
                    f_U_anterior=costo_lineas_anterior+penalizacion_anterior; 
                     if (~isempty(Edges_anterior{1,1}))                   
                        Edges_anterior(1,2:end)={costo_lineas_anterior penalizacion_anterior solucion_compe_anterior1};
                        MapEdgesLocal_2(Edges_anterior{1,1})=[Edges_anterior{1,4}];                    
                        MapEdgesLocal(Edges_anterior{1,1}) = [Edges_anterior{1,2} Edges_anterior{1,3} ];                    
                     end 
        %% Aplicar selección (DE)              
                    if f_U_anterior<f_costo_X_anterior1
                        f_costo_X_anterior1=f_U_anterior;
                        poblacion_anterior1=particle_mod_anterior;
                        Comp_Indu_Capa_etap_anterior1=solucion_compe_anterior1;
                    else
                        Comp_Indu_Capa_etap_anterior1=Comp_Indu_Capa_etap1_anterior;                
                    end
                end

    %% Actualizar el numero de evaluaciones de la funcion objetivo                     
                Eva_BL_anterior=Eva_BL_anterior+Eva;    
            end
        end
    end         
end