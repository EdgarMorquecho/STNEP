function [Gen_Pot_Act_Tot,Costo_Total_Generacion,Comp_Shunt_Nodo,...
        Barra_con_shunt,solucion_compe]=Calculo_Comp_Gen_Perd(x_posicion, ...
        sistema,xopt,~,~,...
        ~,N_exp,incluir_contingencias,~)
  
[Prototi_Contin]=Nodos_Contin(sistema,incluir_contingencias); % Cargar las líneas a las que se aplicar contingencias 
%% Determinar el numero y la posicion de la contingencia    
indx_contingencia = find(Prototi_Contin~=0); % Identificar posicion donde se adicionan lineas de transmision para el individuo
Numero_contingencias=numel(indx_contingencia);
mpc=loadcase(sistema); 
Comp_Indu_Capa=zeros(size(mpc.bus,1),1); %variable para guardar la compensacion para el siguietne a?o
Costo_shunt_total_barra=zeros(size(mpc.bus,1),1); %variable para guardar el costo de la compensacion para el siguietne a?o

%% Determinar la dimension del sistema    
Gen_real_y_ficticio=mpc.gen(:,1); % Lectura de barras con generadores (reales y ficticios)
Size_Gen_real_y_ficticio=size(Gen_real_y_ficticio); %Determinar el numero de generadores reales y ficticios
[Posicion_Nodos_unicos,Nodos_unicos]=unique(Gen_real_y_ficticio,'rows');%Determina barras no repetidas
nodos_repetidos=setdiff(1:Size_Gen_real_y_ficticio(1),Nodos_unicos);%Identifia posicion de los nodos no repetidos
valores_repetidos=Gen_real_y_ficticio(nodos_repetidos, 1); %Identifica posicion de los nodos repetidos  
ind_flag=ismember(Gen_real_y_ficticio, valores_repetidos); %Coloca cero en la posicion de nodos no repetidos
Total_Compensa=zeros(Size_Gen_real_y_ficticio(1),Numero_contingencias+1); % matriz para almacenar la compensacion shunt (condiciones normales y contingencias)
barra_identificada=1;
 for j=1:Size_Gen_real_y_ficticio(1) 
    if ind_flag(j,:)>0 && barra_identificada==1
        Total_Compensa(j,:)=-inf; % colocar -inf para luego calcular la compensación máxima shunt
        barra_identificada=barra_identificada+1;
    elseif ind_flag(j,:)>0 && barra_identificada==2
        Total_Compensa(j,:)=inf; % colocar inf para luego calcular la compensación mínima shunt
        barra_identificada=1;
    else 
        barra_identificada=1;
    end
 end 
%% Inicializa las penalizaciones para el caso base y contingencia
   penalizacion=0;%% Inicializa la penalizacion total
     nc=1; % inicializar contador de contingencias  (uso de gen, artificiales)     
%% Lazo para analizar el caso base y contingencias
for nl=1:1+Numero_contingencias
%% Inicializa las penalizaciones para el caso base y contingencia
        penalizacion1=0;
        penalizacion2=0;
            linea_i_j_total= mpc.branch;     
%% Calculo de las parametros del sistema (resistencia, reactancia, susceptancia) 
        if nl==1 % Caso base
%% Calculo de las parametros del sistema (resistencia, reactancia, susceptancia)
            for j=1:size(x_posicion,1)
                if x_posicion(j,1)>0
                    linea_j=mpc.ne_branch(j,1:end-1);
                    linea_j(1,11)=1;
                    linea_j_adicionadas=repmat(linea_j,x_posicion(j,1),1);
                    linea_i_j_total=[linea_i_j_total;linea_j_adicionadas];
                end
            end
            mpc.branch=linea_i_j_total;   
        else
            for j=1:size(x_posicion,1)
                if x_posicion(j,1)>0
                    linea_j=mpc.ne_branch(j,1:end-1);
                    linea_j(1,11)=1;
                    linea_j_adicionadas=repmat(linea_j,x_posicion(j,1),1);
                    linea_i_j_total=[linea_i_j_total;linea_j_adicionadas];
                end
            end
            mpc.branch=linea_i_j_total;    
            indx_contingencia1=ismember(mpc.branch(:,1:2),mpc.ne_branch(indx_contingencia(nl-1),1:2),'rows');
            f=find(indx_contingencia1>0);
            mpc.branch(f(1,1),:) = [];                       
        end    

%% CALCULAR FLUJO OPTIMO DE POTENCIA CON MATPOWER version 6
        opt=mpoption('OPF_ALG', 540, 'OUT_ALL', 0);            
        [results, success]=runopf(mpc, opt); 
        if nl==1
%% CALCULO DE COSTO DE PERDIDAS          
    %% CALCULO DE PERDIDAS (para cada estado)

    %% CALCULO DE COSTO DE OPERACION (para cada estado)
            Genera=results.gencost(1:Size_Gen_real_y_ficticio(1),5:7);% Lectura del costo de generadores ficticios de potencia activa             
            Genera((ind_flag(:,1)==1),:)=[]; % Eliminar costo de generadores ficticios  
            Gen_barra_actual=results.gen(:,2); % Produccion de MW de todos los genradores (generadores y generadores ficticios)
            Gen_barra_actual((ind_flag(:,1)==1),:)=[]; % Eliminar generadores ficticios 
            Gen_Pot_Act_Tot=sum(Gen_barra_actual);% Generacion de potencia activa total actual            
            %% Calculo del costo de operacion
            Costo_por_generacion=Genera; % Costo de generadores
            Ecuacion=3;
            Costo_generacion=zeros(size(Gen_barra_actual,1),3);
             for j=1:Ecuacion
                 Ecuacion=Ecuacion-1;
                 Costo_generacion(:,j)=Costo_por_generacion(:,j).*Gen_barra_actual.^(Ecuacion);        
             end
            [Fp]=Factor_plant(sistema); % Extraer el factor de planta de cada generador
            Costo_generacion_actual=Fp.*Costo_generacion(:,2); % Costo de generación de cada generador
            Costo_Generacion_actual=sum(Costo_generacion_actual)+sum(Costo_generacion(:,3));    % suma del costo total de cada generador real
            Costo_Total_Generacion=sum(Costo_Generacion_actual); % costo total de generación (costo al valor presente)    
        end 

%% Calculo de la produccion activa de los generadores ficticios
        Costo_Generacion_ficticia=results.gencost(1:Size_Gen_real_y_ficticio(1),6);
        Costo_Generacion_ficticia((ind_flag(:,1)~=1),:)=[]; % Eliminar costo de generadores no ficticios
        G_ficticia_barra=results.gen(:,2); % Produccion de MW de todos los generadores (generadores y generadores ficticios)
        G_ficticia_barra((ind_flag(:,1)~=1),:)=[]; % Eliminar generadores no ficticios
        Costo_Generacion_ficticia_producida= sum(Costo_Generacion_ficticia.*G_ficticia_barra);               
%% CALCULO DE LA COMPENSACIÓN SHUNT
     Gen_shunt_reactiva_actual=results.gen(:,3);
         for j=1:Size_Gen_real_y_ficticio(1)
             if ind_flag(j,1)==0
                 Gen_shunt_reactiva_actual(j,:)=0;
             end
         end 
     Total_Compensa(:,nl)=Gen_shunt_reactiva_actual; % almacenar compensación shunt para la condición base      
%% VERIFICAR LA FACTIBILIDAD DE CADA TOPOLOGIA
 %% Si usa generadores ficticios
        if ((abs(Costo_Generacion_ficticia_producida))>=0.001)&&(success==1)% topologia presenta corte de carga
            penalizacion1=10E7*round((abs(Costo_Generacion_ficticia_producida)));
            nc=nc+1; % contador de topologías que presentan corte de carga        
        end
 %% Si la topologia no converge      
        if success==0 % topologia no es factible 
            penalizacion2=10E15;
            penalizacion=penalizacion+penalizacion1+penalizacion2; % suma de penalizaciones
            nc=1; % inicializar contador de contingencias
            break
        else
            penalizacion=penalizacion+penalizacion1+penalizacion2; % suma de penalizaciones            
        end
end
%% DETERMINAR LA MÁXIMA COMPENSACIÓN (CONDICIÓN NORMAL Y BAJO CONTINGENCIAS)
C_max_Cont=zeros(Size_Gen_real_y_ficticio(1),1); %Matriz para almacenar la compensacion shunt maxima entre el caso base y bajo contingencias
barra_identificada=1;
         for j=1:Size_Gen_real_y_ficticio(1) 
            if ind_flag(j,:)>0 && barra_identificada==1
                C_max_Cont(j,:)=max(Total_Compensa(j,:),[],2); % halla la maxima compensacion en cada nodo (compensacion capacitiva)
                barra_identificada=barra_identificada+1;
            elseif ind_flag(j,:)>0 && barra_identificada==2
                C_max_Cont(j,:)=min(Total_Compensa(j,:),[],2); % halla la compensacion shunt minima en cada nodo (compensacion indictiva)
                barra_identificada=1;
            else 
                barra_identificada=1;
            end
         end
Comp_Shunt_Nodo=C_max_Cont;
Barra_con_shunt=ind_flag;
%% GUARDAR LA COMPENSACIÓN SHUNT ENCONTRADA PARA ADICIONAR AL SIGUIENTE ESTADO
    costo_compensacion_actual=results.gencost(Size_Gen_real_y_ficticio(1)+1:2*Size_Gen_real_y_ficticio(1),6); %separar los costo de generacion de MW de MVAr
    for j=1:size(results.bus,1)
    Comp_ficticia=find(Gen_real_y_ficticio==Posicion_Nodos_unicos(j));
        if size(Comp_ficticia,1)==1
            Comp_Indu_Capa((j),1)=0; 
        elseif size(Comp_ficticia,1)==2
            Comp_Indu_Capa((j),1)=(C_max_Cont(Comp_ficticia(1),1)+C_max_Cont(Comp_ficticia(2),1));
            if Comp_Indu_Capa((j),1)>=0
                Costo_shunt_total_barra((j),1)=Comp_Indu_Capa((j),1)*costo_compensacion_actual(Comp_ficticia(1),1);  
            else
                Costo_shunt_total_barra((j),1)=Comp_Indu_Capa((j),1)*costo_compensacion_actual(Comp_ficticia(2),1);    
            end       
        end
    end                                             
    solucion_compe=Comp_Indu_Capa(:,1); % guarda la compensación máxima de la topología actual (condición normal y contingencias)
   end
