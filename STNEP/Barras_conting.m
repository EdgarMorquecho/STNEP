function [Contingencias_entre_barras]=Barras_conting(x_posicion,sistema,Contin_Segun_Dimensi,Contin_Segun_Dimensi_Final,incluir_contingencias)
%    global sistema
%    global mpc
 [resistencia_total, reactancia_total, susceptancia_total, limite_flujo, modelo]=Topologia_Red(x_posicion, sistema);
[LineasIniciales, ~]=TestSystems(sistema); % recoje los costos y la 
% % topologia inicial
numero_ramas=numel(LineasIniciales); 
% 
mpc=loadcase(sistema);
define_constants;

% %%% Addecuacion del status de las lineas
mpc.branch(:,4)=reactancia_total;
mpc.branch(:,3)=resistencia_total;
mpc.branch(:,5)=susceptancia_total;
mpc.branch(:,6)= limite_flujo;

Contingencias_entre_barras=char(numero_ramas,1);
Se_produjo_contingencia=char('C');
No_produjo_contingencia=char('N');
    for j=1:numero_ramas
        Contingencias_entre_barras(j)=No_produjo_contingencia;
    end
status=LineasIniciales'|x_posicion;
mpc.branch(:,11)=status;
inicial_mpc=mpc;
%% Inicializa las penalizaciones
    penalizacion1=0;
    penalizacion2=0;
%% Inicializa las penalizaciones
%% Chequear si no es nodo aislado
switch modelo
        case {'AC'}

    opt=mpoption('OPF_ALG', 540, 'OUT_ALL', 0);
    [results, success]=runopf(mpc, opt);
    
%% CALCULO DE COSTO DE OPERACION (para cada estado) 
    shunt=results.gen(:,1); % Identificar entre generador y generador ficticio
    Ada=size(shunt);
    nodo_repetido=shunt;
    [~,compensacion_unicos]=unique(nodo_repetido,'rows');
    nodos_repetidos=setdiff(1:Ada(1),compensacion_unicos) ;
    valores_repetidos=shunt(nodos_repetidos, 1);   
    ind_flag=ismember(shunt, valores_repetidos);
         Genera=results.gencost(1:Ada(1),5:7);
         for fg=1:Ada(1)
             if ind_flag(fg,1)==1
                 shunt(fg)=-1; 
             end
         end
              
     Genera((shunt(:,1)==-1),:)=[];  % Eliminar costo de generadores ficticios
     G_barra=results.gen(:,2);       % Produccion de MW de todos los genradores (genradores y genradores ficticios)
     G_barra((shunt(:,1)==-1),:)=[]; % Eliminar generadores ficticios
     generacion_barra=G_barra;       % Produccion  total de MW de cada generador
     Costo_por_generacion=Genera;    % Costo  de generadores
     Ecuacion=3;
     Costo_generacion=zeros(size(generacion_barra,1),3);
     for hj=1:Ecuacion
         Ecuacion=Ecuacion-1;
         Costo_generacion(:,hj)=Costo_por_generacion(:,hj).*generacion_barra.^(Ecuacion);        
     end
     C_T=sum(Costo_generacion);
     C_T_Generacion=sum(C_T); % separar costo de operacion de penalizacion 
     
  %% CALCULO DEL COSTO DE LA COMPENSACION REACTIVA     
    costo_100_reactiva=results.gen(:,3);
    reactivacompen_consto=results.gencost(Ada(1)+1:2*Ada(1),6); %separar los costo de generacion de MW de MVAr
         for ghg=1:Ada
             if ind_flag(ghg,1)==0
                 costo_100_reactiva(ghg,:)=0;
             end
         end

     reactiva_cost=costo_100_reactiva.*reactivacompen_consto; % calcular el costo de cada compensador
     reactiva_cost_total=sum(reactiva_cost); % calcular la compensacion shunt total para cad estado     
     
     
%% VERIFICAR LA FACTIVILIDAD DE CADA TOPOLOGIA     
   C_sin_generacion=results.f-C_T_Generacion-reactiva_cost_total; % verifiar si existe generacion  ficticia (penalizacion)
    if (round(abs(C_sin_generacion))>=0)&&(success==1)% topologia presenta corte de carga
        penalizacion1=1000*round((abs(C_sin_generacion))); % penalizar por corte de carga
    end
        
    if success==0 % topologia no es factible 
        penalizacion2=abs(C_sin_generacion)+10E9; % penalizar por solucion no factible
    end
    penalizacion=penalizacion1+penalizacion2; % suma de penalizaciones    
                  
    
 if incluir_contingencias==1    
 %CONTINGENCIA        
 if penalizacion==0%topologia actual no presenta corte de carga
   if Contin_Segun_Dimensi==0 
    for barra=1:numero_ramas
       if LineasIniciales(barra)>0
           Contingencias_entre_barras(barra,1)=Se_produjo_contingencia;
       else
           Contingencias_entre_barras(barra,1)=Contingencias_entre_barras(barra,1);
       end
        
    end

     else
          [Prototi_Contin]=Nodos_Contin(sistema,incluir_contingencias);
           for k=1:numero_ramas
               if Prototi_Contin(k)>0
               if x_posicion(k)>0
                  Contingencias_entre_barras(k,1)=Se_produjo_contingencia;                                
                else
                  Contingencias_entre_barras(k,1)=Contingencias_entre_barras(k,1);
               end
               end

            end
   end
 end
 end
end