function [poblacion_anterior1,f_costo_X_anterior1,Comp_Indu_Capa_etap_anterior1,Eva_BL_anterior]...
                =Elimina_linea_innecesaria_cont(poblacion_anterior1,Xmin,~,sistema,...
                incluir_contingencias,incluir_perdidas,factor_perdidas,...
                compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
                SS,~,...
                Comp_Indu_Capa_etap_anterior1,f_costo_X_anterior1,...
                solucion_compe_anterior1,...
                ~,fhd,Contin_Segun_Dimensi_Final,MapEdgesLocal_2,MapEdgesLocal,Costos)
            

%%% ELIMINAR LINEAS ADICIONADAS QUE POSIBLEMENTE SON INNECESARIAS
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%% Identificar las lineas adicionadas para el individuo actual en cada etapa 
Eva_BL_anterior=0;
    Cadd=poblacion_anterior1; %Identificar las lineas adicionadas en el individuo para (n_anos=1)                                                                                                     
    indx = find(Cadd~=0); % Identificar posicion donde se adicionan lineas de transmision para el individuo
%% Inicio del proceso iterativo                    
    if ~isempty(indx)  %Elimina lineas en funcion de las perdidas si existen lineas adicionadas en el individuo actual
        indices_perm=randperm(size(indx,1));
%% Iniciar proceso de eliminacion de lineas        
        for indx_j=1:size(indx,1)
                Si_cambio=0;
            for indx_k=1:Cadd(indx(indices_perm(indx_j)))
                if Si_cambio==0
                particle_mod_anterior = poblacion_anterior1;
    %% Eliminar la linea indx_j de lineas adicionadas
                particle_mod_anterior(indx(indices_perm(indx_j)),1)=particle_mod_anterior(indx(indices_perm(indx_j)),1)-(Cadd(indx(indices_perm(indx_j)))-(indx_k-1));%Eliminar una linea 
                Edges_anterior = cell(1,4);
                Eva=0;  
                Comp_Indu_Capa_etap1_anterior=Comp_Indu_Capa_etap_anterior1;
%% Determinar si el individuo actual presenta mejorvalor por adicion de lineas que el individuo anterior      
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
                        Si_cambio=1;
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

end