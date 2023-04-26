function [costo_total, penalizacion,solucion_compe]=evaluacion_particulas_aux(poblacion,sistema,Qp,Incremen_demanda_final_React,...
    Incremen_demanda_final_Act,Incremen_generacion_final, n_anos,compensacion_anual_optim,Comp_Indu_Capa_inicial,...
    Interes_etapa,Comp_Indu_Capa_etap,...
       Incr_Demam_React,Incr_Demam_Act,Incremento_Gen_Q_max,Incremento_Gen_Q_min,Incremento_Gen_P_max,Incremento_Gen_P_min)

mpc=loadcase(sistema);
define_constants;
solucion_compe=0;

%% Lee la topologia inicial para el estado n

%Solo para eliminación de lineas
aux=strcmp(Qp, 'eliminacion_lineas');
if aux==1
    Qp=1;
    x_posicion=poblacion(:,Qp);
    [~,Costos]=TestSystems(sistema); % recoje los costos y la topologia inicial del estado n-1
    LineasIniciales=poblacion(:,Qp)';%(:,i-SS/n_anos)';
else
    x_posicion=poblacion(:,Qp);
    %Dinámico
    if Qp>1
        [~,Costos]=TestSystems(sistema); % recoje los costos y la topologia inicial del estado n-1
        %LineasIniciales=Lineas_ano(:,i-SS/n_anos)';
        LineasIniciales=poblacion(:,Qp-1)';%(:,i-SS/n_anos)';
        %Estático
    else
        [LineasIniciales,Costos]=TestSystems(sistema); % recoje los costos y la topologia inicial del sistema
    end
end

numero_ramas=numel(LineasIniciales);
LT_nuevas=x_posicion-LineasIniciales';


%% Calculo de las parametros del sistema (resistencia, reactancia, susceptancia)
[resistencia_total, reactancia_total, susceptancia_total, limite_flujo]=Topologia_Red_aux(poblacion,sistema,Qp);

%  Actualizacion del estado y propiedades de las lineas
mpc.branch(:,4)=reactancia_total;
mpc.branch(:,3)=resistencia_total;
mpc.branch(:,5)=susceptancia_total;
mpc.branch(:,6)= limite_flujo;
status=LineasIniciales'|x_posicion;
mpc.branch(:,11)=status;


%% Inicializa las penalizaciones
penalizacion1=0;
penalizacion2=0;

%% Incremento de la demanda activa y reactiva (incremento porcentual)
%[Incr_Demam_React,Incr_Demam_Act]=incremen_deman(Inc_Demanda_per_cent,sistema,Qp);
Incr_Demam_React=Incremen_demanda_final_React;
Incr_Demam_Act=Incremen_demanda_final_Act;
mpc.bus(:,4)=Incr_Demam_React;
mpc.bus(:,3)=Incr_Demam_Act;


%% Incremento de la generacion activa y reactiva (incremento porcentual)
%[Incremen_generacion_final]=incremen_generacion(Inc_Demanda_per_cent,sistema,Qp);
mpc.gen(:,4)=Incremento_Gen_Q_max(:,Qp);
mpc.gen(:,5)=Incremento_Gen_Q_min(:,Qp);
mpc.gen(:,9)=Incremento_Gen_P_max(:,Qp);
mpc.gen(:,10)=Incremento_Gen_P_min(:,Qp);


%% Inserta la compensacion shunt del estado anterior en el estado actual
if compensacion_anual_optim==1
    
    %inicio de variables para compensacion
    Comp_Indu_Capa=zeros(size(mpc.bus,1),n_anos); %variable para guardar la compensacion para el siguietne a?o
    Costo_shunt_total_barra=zeros(size(mpc.bus,1),1); %variable para guardar el costo de la compensacion para el siguietne a?o
    
    if Qp==1
        mpc.bus(:,6)=Comp_Indu_Capa_inicial; % compensacion shunt para el estado 1
    else
        Com_trf=Comp_Indu_Capa_inicial; % compensacion shunt para el estado diferente de 1
        for ttgb=1:Qp-1
            Com_trf=Com_trf+Comp_Indu_Capa_etap(:,ttgb);% compensacion shunt para el estado n+1
            %Com_trf=Comp_Indu_Capa_etap(:,ttgb);% compensacion shunt para el estado n+1
        end
        mpc.bus(:,6)= Com_trf;
    end
end


%% CALCULAR FLUJO OPTIMO DE POTENCIA CON MATPOWER
mpopt = mpoption('model','AC', 'opf',struct('ac',struct('solver','PDIPM'),'current_balance',0,'v_cartesian',0, 'violation',5e-6,'use_vg',0,'flow_lim','S','ignore_angle_lim',0),'out.all', 0, 'pdipm.feastol',0);     
[results, success]=runopf(mpc, mpopt);

if compensacion_anual_optim==1
    
    %         flujo_sistema_max=results.branch(:,6);
    %         %flujo_sistema=abs([results.branch(:,14) results.branch(:,16)])
    %         flujo_sistema=max(abs([results.branch(:,14) results.branch(:,16)]),[],2);
    %         Indice_flujo1=zeros(size(flujo_sistema_max,1),1);
    %         for j=1:numero_ramas
    %             if  flujo_sistema_max(j,1)==0
    %                 Indice_flujo1(j,1)=0;
    %             else
    %                 Indice_flujo1(j,1)=flujo_sistema(j,1)*100/flujo_sistema_max(j,1);
    %             end
    %         end
    %         Indice_flujo2(:,1)=Indice_flujo1; % almacenar compensación shunt para la condición base
    %
    %% CALCULO DE COSTO DE OPERACION (para cada estado)
    shunt=results.gen(:,1); % Identificar entre generador y generador ficticio
    Ada=size(shunt);
    nodo_repetido=shunt;
    [~,compensacion_unicos]=unique(nodo_repetido,'rows');
    nodos_repetidos=setdiff(1:Ada(1),compensacion_unicos);%Identifia posicion de los nodos no repetidos
    valores_repetidos=shunt(nodos_repetidos, 1); %Identifica posicion de los nodos repetidos
    ind_flag=ismember(shunt, valores_repetidos); %Coloca cero en la posicion de nodos no repetidos
    %Genera=results.gencost(1:Ada(1),5:7);
    %Numero_gen=mpc.gen(1:Ada(1),1);%
    for j=1:Ada(1)
        if ind_flag(j,1)==1
            shunt(j)=-1;
        end
    end
    
    %     Genera((shunt(:,1)==-1),:)=[];  % Eliminar costo de generadores ficticios
    %     Numero_gen((shunt(:,1)==-1),:)=[];%
    %     G_barra=results.gen(:,2);       % Produccion de MW de todos los genradores (genradores y genradores ficticios)
    %     G_barra((shunt(:,1)==-1),:)=[]; % Eliminar generadores ficticios
    %     generacion_barra=G_barra;       % Produccion  total de MW de cada generador
    %     Costo_por_generacion=Genera;    % Costo  de generadores
    %     Ecuacion=3;
    %     Costo_generacion=zeros(size(generacion_barra,1),3);
    %     for j=1:Ecuacion
    %          Ecuacion=Ecuacion-1;
    %          Costo_generacion(:,j)=Costo_por_generacion(:,j).*generacion_barra.^(Ecuacion);
    %      end
    %      C_T=sum(Costo_generacion);
    %      C_T_Generacion=sum(C_T); % separar costo de operacion de penalizacion
    
    
    
    
    %% CALCULO DEL COSTO DE LA COMPENSACION REACTIVA
    costo_100_reactiva=results.gen(:,3);
    reactivacompen_consto=results.gencost(Ada(1)+1:2*Ada(1),6); %separar los costo de generacion de MW de MVAr
    for j=1:Ada
        if ind_flag(j,1)==0
            costo_100_reactiva(j,:)=0;
        end
    end
    
    %reactiva_cost=costo_100_reactiva.*reactivacompen_consto; % calcular el costo de cada compensador
    %reactiva_cost_total=sum(reactiva_cost); % calcular la compensacion shunt total para cad estado
    
    
    % %      %% VERIFICAR LA FACTIVILIDAD DE CADA TOPOLOGIA
    % %    C_sin_generacion=results.f-C_T_Generacion-reactiva_cost_total; % verifiar si existe generacion  ficticia (penalizacion)
    % %     if (round(abs(C_sin_generacion))>=0)&&(success==1)% topologia presenta corte de carga
    % %         penalizacion1=1000*round((abs(C_sin_generacion))); % penalizar por corte de carga
    % %     end
    % %
    % %     if success==0 % topologia no es factible
    % %         penalizacion2=abs(C_sin_generacion)+10E9; % penalizar por solucion no factible
    % %     end
    % %     penalizacion=penalizacion1+penalizacion2; % suma de penalizaciones
    
    %% CALCULO DEL COSTO DE COMPENSACIÓN
    Ale=results.gen(:,3);
    Ale2=results.gen(:,1);
    for j=1:Ada(1)
        if ind_flag(j,1)==0
            Ale(j,:)=0;
            Ale2(j,:)=0;
        end
    end
    Neta=size(mpc.gen);
    Total_Compensa=zeros(Neta(1),numero_ramas+1); % matriz para almacenar la compensacion shunt (condiciones normales y contingencias)
    barra_identificada=1;
    for j=1:Ada(1)
        if Ale2(j,:)>0 && barra_identificada==1
            Total_Compensa(j,:)=-inf; % colocar -inf para luego calcular la compensación máxima shunt
            barra_identificada=barra_identificada+1;
        elseif Ale2(j,:)>0 && barra_identificada==2
            Total_Compensa(j,:)=inf; % colocar inf para luego calcular la compensación mínima shunt
            barra_identificada=1;
        else
            barra_identificada=1;
        end
    end
    
    Total_Compensa(:,1)=Ale; % almacenar compensación shunt para la condición base
    
    %% DETERMINAR LA MÁXIMA COMPENSACIÓN (CONDICIÓN NORMAL Y BAJO CONTINGENCIAS)
    C_max_Cont=zeros(Ada(1),1);
    barra_identificada=1;
    for j=1:Ada(1)
        if Ale2(j,:)>0 && barra_identificada==1
            C_max_Cont(j,:)=max(Total_Compensa(j,:),[],2); % halla la maxima compensacion en cada nodo (compensacion capacitiva)
            barra_identificada=barra_identificada+1;
        elseif Ale2(j,:)>0 && barra_identificada==2
            C_max_Cont(j,:)=min(Total_Compensa(j,:),[],2); % halla la compensacion shunt minima en cada nodo (compensacion indictiva)
            barra_identificada=1;
        else
            barra_identificada=1;
        end
    end
    
    
    %% GUARDAR LA COMPENSACIÓN SHUNT ENCONTRADA PARA ADICIONAR AL SIGUIENTE ESTADO
    if Qp>0
        sn=0;
        Am=0;
        mb=1;
        for ad=1:Ada(1)
            if shunt(ad)~=-1
                sn=1+sn;
                Comp_Indu_Capa(mb,Qp)=0;
                mb=1+mb;
            elseif shunt(ad)==-1
                Am=Am+1;
                if Am<=1
                    sn=2+sn;
                    Comp_Indu_Capa(mb,Qp)=(C_max_Cont(sn-1,1)+C_max_Cont(sn,1));
                    if Comp_Indu_Capa(mb,Qp)>=0
                        Costo_shunt_total_barra(mb,1)=Comp_Indu_Capa(mb,Qp)*reactivacompen_consto(sn-1,1)/(1+Interes_etapa)^(Qp-1);
                    else
                        Costo_shunt_total_barra(mb,1)=Comp_Indu_Capa(mb,Qp)*reactivacompen_consto(sn,1)/(1+Interes_etapa)^(Qp-1);
                    end
                    mb=1+mb;
                else
                    Am=0;
                end
            end
        end
        solucion_compe=Comp_Indu_Capa(:,Qp); % guarda la compensación máxima de la topología actual  (condición normal y contingencias)
    end
    
    
    C_Shunt_maxima_Costo=sum(Costo_shunt_total_barra); %costo de compensación shunt total ( condición normal y bajo contingencias)
    
    
else
    C_Shunt_maxima_Costo=0;
    
end












%% PENALIZACION GEN ARTIFICIALES

indices_gen=mpc.gen(:,1);
%ind_gen_reales=indices_gen;
aux1=ones(length(indices_gen),1);
contador=0;
for aux=1:length(indices_gen)
    if length(find(aux1.*indices_gen(aux)==indices_gen))>1
        contador=contador+1;
        %ibus_gen_artif(contador)=indices_gen(aux);
        i_gen_artif(contador)=aux;
    end
end
%ibus_gen_reales=setdiff(indices_gen,ibus_gen_artif');
i_gen_reales=setdiff(indices_gen,i_gen_artif');
i_gen_artif=unique(i_gen_artif(find(i_gen_artif)))';

gen_artificial_total=sum(round(results.gen(i_gen_artif,2)));

%--Penalizacion de GEN. ARTIFICIALES-------------------
%if (round(abs(results.f))>=0)&&(success==1)

%             if (success==1)&&(round(abs(gen_artificial_total))==0)
%                 1+1;
%             end

if (success==1)&&(round(abs(gen_artificial_total))>0)
    penalizacion1=round(abs(gen_artificial_total))*10E7;    %Gen. Ficticios    results.f=valor de la F.O corresp. a $Gen.Artif. que ayudan a Converger y OPTIMIZAR
end
%------------------------------------------------
%--penalizacion NO convergencia------------------
if success==0
    penalizacion2=abs(gen_artificial_total)+10E17;                 %Penalizacion de NO converegencia
end
%------------------------------------------------
penalizacion=penalizacion1+penalizacion2;


%costo_total=sum(Costos'.*LT_nuevas);%+numCapacitores*0;%+penalizacion;
costo_total=sum((Costos/(1+Interes_etapa)^(Qp-1))'.*LT_nuevas);%+numCapacitores*0;%+penalizacion

%angulo_buses=results.bus(:,9);      %para creacion sol.ini (criterio min esfuerzo)

penalizacion=round(penalizacion);%+C_Shunt_maxima_Costo;
%penalizacion=round(penalizacion)+C_T_Generacion11+C_Shunt_maxima_Costo+Costo_perdidas/(1+Interes_etapa)^(Qp-1);


clearvars -except costo_total penalizacion Eva solucion_compe i x_posicion...
    sistema Qp Comp_Indu_Capa_etap Contin_Segun_Dimensi_Final...
    incluir_contingencias incluir_perdidas Inc_Demanda_per_cent factor_perdidas...
    compensacion_anual_optim landa Lineas_ano Comp_Indu_Capa_inicial...
    Comp_Indu_Capa Interes_etapa SS n_anos angulo_buses Reactancia_Usada...
    Resistencia_Usada Indice_flujo losses_linea Costo_despacho_lineas_min Costo_despacho_lineas_max

end
