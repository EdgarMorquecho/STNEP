function varargout = TNEP_Dinam(varargin)
% TNEP_DINAM MATLAB code for TNEP_Dinam.fig
%      TNEP_DINAM, by itself, creates a new TNEP_DINAM or raises the existing
%      singleton*.
%
%      H = TNEP_DINAM returns the handle to a new TNEP_DINAM or the handle to
%      the existing singleton*.
%
%      TNEP_DINAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TNEP_DINAM.M with the given input arguments.
%
%      TNEP_DINAM('Property','Value',...) creates a new TNEP_DINAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TNEP_Dinam_OpeningFcn gets called. An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TNEP_Dinam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TNEP_Dinam

% Last Modified by GUIDE v2.5 15-Apr-2021 09:25:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TNEP_Dinam_OpeningFcn, ...
                   'gui_OutputFcn',  @TNEP_Dinam_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   %% PARAMETROS DEL PET
   %% Llama al tipo de sistema
   delete(figure(1))
   contents = cellstr(get(handles.popupmenu1,'String'));%Nombre del Sistema ingresado en la interfaz
   sistema=contents{get(handles.popupmenu1,'Value')}; 
   
   %% Poblacion Inicial   
   Opcion_poblacion=cellstr(get(handles.popupmenu6,'String'));% Nombre de la poblacion inicial ingresada en la interfaz
   opcion=Opcion_poblacion{get(handles.popupmenu6,'Value')};     

   %% Maximo numero de lineas por derecho de trasmision
   Xmax_ingresado=get(handles.edit9,'String'); %Maximo numero de lineas por derecho de trasmision ingresada en la interfaz
   Xmax= str2double(Xmax_ingresado);
   set(handles.edit9,'string',Xmax); 
   
   
   %% Tamano de poblacion
   numero_individuos=get(handles.edit3,'String'); %Tamano de poblacion ingresada en la interfaz
   Particle_Number= str2double(numero_individuos); %Multiplicar el tamano de la poblacion por el numero de estados dinamicos totales
   Particle_Number_original= str2double(numero_individuos);
   set(handles.edit3,'string',Particle_Number_original); 
   
   %% Valor optimo de la funcion objetivo en caso de que se lo conozca para sistemas de prueba 0 por defecto
   Valor_referencial=get(handles.edit4,'String'); %Valor de referencia (minimo) si se conoce caso contrario cero por defecto
   Valor_min_referencia= str2double(Valor_referencial);
   set(handles.edit4,'string',Valor_min_referencia);   
   
   %% Maximo numero de iteraciones o generaciones (criterio de parada)
   Maxima_iteracion=get(handles.edit5,'String') ;
   Max_Gen= str2double(Maxima_iteracion);
   set(handles.edit5,'string',Max_Gen);    
   
   %% Numero de veces que se quiere repetir el test
   Ensayos_ingresado=get(handles.edit6,'String');
   nexp= str2double(Ensayos_ingresado);
   set(handles.edit6,'string',nexp); 
   
   %% Numero de iteraciones en las que se repite el valor encontrado sin cambiar
   Valor_min_repetido=get(handles.edit12,'String');
   Numero_iter_sin_cambiar= str2double(Valor_min_repetido);
   set(handles.edit12,'string',Numero_iter_sin_cambiar);    

   %% Incluir contingencias
   %% Ingresar datos para aplicar o no contingencias
   Contingencias_ingresado=cellstr(get(handles.popupmenu2,'String'));
   security=Contingencias_ingresado{get(handles.popupmenu2,'Value')};   
   valor_contingencia=strcmp(security, 'Yes');
   if valor_contingencia==1
     incluir_contingencias=1; % Si es 1 realizar contingencias
   else
      incluir_contingencias=0; %Si es cero no realizar contingencias;
   end
%% Ingresar datos para aplicar contingencias al sistema base o un vector de contingencias
   Vector_contingencias_ingresado=cellstr(get(handles.popupmenu8,'String'));
   Aplicar_contingencia=Vector_contingencias_ingresado{get(handles.popupmenu8,'Value')};   
   Contingencia_vector_o_initial_topologia=strcmp(Aplicar_contingencia, 'Contingency_vector');
   if Contingencia_vector_o_initial_topologia==1
      Contin_Segun_Dimensi_Final=1; % Si es 1 realizar contingencias usando el vector de contingencias
   else
      Contin_Segun_Dimensi_Final=0; %Si es cero realiza contingneicas al sistema base;
   end
   
%% Parametros relacionados a las perdidas  
   %Incluir perdidas
   losses = cellstr(get(handles.popupmenu3,'String'));
   losses1=losses{get(handles.popupmenu3,'Value')};   
   valor_perdidas=strcmp(losses1, 'Yes');
   if valor_perdidas==1
     incluir_perdidas=1; % Si es =1 incluir costo de perdidas
   else
      incluir_perdidas=0; % Si es =1 no incluir costo de perdidas
   end    
   
   %% Valor de lamda (perdidas)
   Costo_lamda=get(handles.edit8,'String') ;
   landa= str2double(Costo_lamda);
   set(handles.edit8,'string',landa);
   
   %% Factor de perdidas (perdidas)
   Perd_factor=get(handles.edit7,'String') ;
   factor_perdidas= str2double(Perd_factor);
   set(handles.edit7,'string',factor_perdidas);     
        
   %% Compensacion dinamica
   Compensacio_dimanica = cellstr(get(handles.popupmenu5,'String'));
   Compensacio_dimanica1=Compensacio_dimanica{get(handles.popupmenu5,'Value')};   
   Compensacio_dimanica2=strcmp(Compensacio_dimanica1, 'Yes');
   if Compensacio_dimanica2==1
     compensacion_anual_optim=1; % incluir compensacion shunt para cada etapa 
   else
      compensacion_anual_optim=0; % no incluir compensacion shunt en cada etapa 
   end
        
   
   %% Maximo numero de evaluaciones de la funcion objetivo
   Max_FES1=get(handles.edit13,'String') ;
   Max_FES= str2double(Max_FES1);
   set(handles.edit13,'string',Max_FES);  
   
%% PARAMETROS ADICIONALES 
   mpc=loadcase(sistema); % cargar el sistema inicial
   dim=size(mpc.ne_branch,1); % Dimension del problema
   
%% Lineas maxima por derecho trasmision   
   Xmax = repmat(Xmax, dim, Particle_Number); % Numero maximo de lineas que se puede adicionar al problema


%% Lineas minimo por derecho trasmision         
   Xmin=zeros(dim,1);
%% evaluacion_particulas: Es la funcion usada para resolverl AC power flow  
   fhd=@evaluacion_particulas; 
   
%%  Sentencia para contar inicializar el tiempo que demora en resolver el problema de planeamiento 
   tic       
   contents1 = cellstr(get(handles.popupmenu4,'String'));
   metaheuristica_used=contents1{get(handles.popupmenu4,'Value')};
     
%% Llamada a la rutina de optimizacion      
    switch metaheuristica_used % Seleccion de la metaheuristica usada para el proceso de optimizacion
     case {'DE_PBILc'}  
        Metaheuristica_Usada='DE_PBILc';
       [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]=DE_PBILc(...
       Max_Gen,Max_FES,Particle_Number,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
       incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
       landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion);        
     case {'IDE_PBILc'}  
        Metaheuristica_Usada='IDE_PBILc';
       [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]=IDE_PBILc(...
       Max_Gen,Max_FES,Particle_Number,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
       incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
       landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion);   
  
     case {'LPSO'}  
        Metaheuristica_Usada='LPSO';
       [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]=LPSO(...
       Max_Gen,Max_FES,Particle_Number,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
       incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
       landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion);     

     case {'TLBO'}  
        Metaheuristica_Usada='TLBO';
       [xopt,fxopt,fevalcount, FOPT, RunStats,C_FINAL_CADA_ETAPA,Comp_Indu_Capa_inicial,iteraciones_cada]=TLBO(...
       Max_Gen,Max_FES,Particle_Number,dim,Xmin,Xmax,Valor_min_referencia, nexp, fhd, sistema,...
       incluir_contingencias,incluir_perdidas,Contin_Segun_Dimensi_Final,factor_perdidas,compensacion_anual_optim,...
       landa,Numero_iter_sin_cambiar,Metaheuristica_Usada,Particle_Number_original,opcion);    
    end

Tiempo_total_de_simlacion=toc; % Tiempo que se tarda para resolver el problema de planeamiento

Comp_final_guardado_etapa=sum(C_FINAL_CADA_ETAPA);
%% CALCULO DEL NUMERO DE LINEAS ADICIONADAS EN CADA PERIODO DINAMICO   
Lineas_adicionadas_cada_estado1=xopt;
   Costo_Lineas_cada_estado1=zeros(dim,nexp);
   mpc=loadcase(sistema);
    for i=1:nexp
        Lineas_adicionadas_cada_estado1(:,i)=xopt(:,i);
        Costo_Lineas_cada_estado1(:,i)=mpc.ne_branch(:,14).*Lineas_adicionadas_cada_estado1(:,i); 
    end
    %----------------------------------
mpc=loadcase(sistema); 
Barras_branch_col=mpc.ne_branch(:,1:2);
Lineas_adicionadas_cada_estado=[Barras_branch_col Lineas_adicionadas_cada_estado1];
Costo_Lineas_cada_estado=[Barras_branch_col Costo_Lineas_cada_estado1];
Barra_generacion_1=mpc.bus(:,1);

%-----------------------------------

%% CALCULO DEL COSTO DE LINEAS ADICIONADAS EN CADA PERIODO DINAMICO
Costo_total_lineas_cada_estado=sum(Costo_Lineas_cada_estado1); %Calculo del costo total de lineas para cada periodo
Costo_Total_lineas_todos_estados=zeros(1,nexp);
    for i=1:nexp
        Costo_Total_lineas_todos_estados(i)=sum(Costo_total_lineas_cada_estado(i:i)); % calculo del costo por adicion para todos los estados dinamicos
    end
%% CALCULO DEL COSTO DE OPERACION PARA CADA ESTADO DINAMICO
Perdidas_Lineas_cada_estado1=zeros(dim,nexp);% variable  para guardar los resultados de perdidas en el sistema
status_linea=zeros(dim,nexp);% variable para ver que lineas fueron adcionadas al final del proceso de optimizacion
Costo_Gen_Pot_Activa_etapa=zeros(1,nexp);
Costo_Gen_Pot_Activa_todas_etapas=zeros(1,nexp);
Gen_Pot_Activa_etapa=zeros(1,nexp);
Gen_Pot_Activa_todas_etapas=zeros(1,nexp);
     mpc=loadcase(sistema);  % carga algunas variables necesarias del sistema de prueba bajo estudio esta es una función de matpower              
     Compensacion_actual_verificada=zeros(size(mpc.bus,1),nexp); %variable para guardar la compensacion para el siguietne a?o
     Comp_shunt_en_cada_barra1=zeros(size(mpc.bus,1),nexp); 
     Costo_individual_compe=zeros(size(mpc.gen,1),nexp);
     Barra_con_shunt=zeros(size(mpc.gen,1),nexp);
for N_exp=1:nexp
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
Perdidas_Lineas_cada_estado1=0;
         
%% CALCULO DE LAS PERDIDAS TOTALES Y EL COSTO PARA CADA NUMERO DE EXPERIMENTO (nexp)
Perdidas_total_cada_estado=sum(Perdidas_Lineas_cada_estado1); % Sumar las perdidas de cada derecho de trasmision (Perdidas totales)
Perdidas_Lineas_cada_estado=0;
Perdidas_total_todos_estados=zeros(1,nexp);
Cost_Perdidas_total_todos_estados=zeros(1,nexp);
Costo_Perdida_cada_estado=zeros(1,nexp);
for i=1:nexp
    Perdidas_total_todos_estados(i)=0;
      Cost_Perdidas_total_todos_estados(i)=0;
end

%% CALCULO DE LA COMPENSACION DE REACTIVOS PARA UNA TOPOLOGIA ESPECIFICA EN CADA ESTADO 
     mpc=loadcase(sistema);  % carga algunas variables necesarias del sistema de prueba bajo estudio esta es una función de matpower              
     Comp_Indu_Capa_barra=zeros(size(mpc.bus,1),nexp); %variable para guardar la compensacion para el siguietne a?o
     Comp_Inductor=zeros(size(mpc.bus,1),nexp); %variable para guardar la compensacion para el siguietne a?o
     Comp_Capacitor=zeros(size(mpc.bus,1),nexp); %variable para guardar la compensacion para el siguietne a?o
     Costo_shunt_total_barra_result=zeros(size(mpc.bus,1),nexp); %variable para guardar el costo de la compensacion para el siguietne a?o
    [Posicion_Nodos_unicos,~]=unique(mpc.gen(:,1),'rows');%Determina barras no repetidas
     for i=1:nexp
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
  Comp_shunt_costo_total=zeros(1,nexp); 
  Compensacion_shunt_total_Stage=zeros(1,nexp); 
  Compensacion_shunt_total=sum(Comp_shunt_en_cada_barra1);% Compensacion shunt total para cada estado  
  for i=1:nexp
      Comp_shunt_costo_total(1,i)=sum(Comp_shunt_costo_etapa(1,i:i));
      Compensacion_shunt_total_Stage(i)=sum(Compensacion_shunt_total(i:i));      
  end
  Comp_shunt_en_cada_barra=[Barra_generacion_1 Comp_shunt_en_cada_barra1];
   
%% VARIABLES NECESARIAS PARA GUARDAR EN UN DOCUMENTO .MAT  
    Particle_name=num2str(Particle_Number_original);  % Transforma el nUmero de individuos en caracteres
    PERDAS=num2str(incluir_perdidas);  % Transforma la variable perdidas en caracter
    SI_continge=num2str(incluir_contingencias);  % Transforma la variable contingenica en 
    E_dinamicas=num2str(1);%Nmero de estados dinamicos
    
%% HORA EN LA QUE FINALIZA LA SIMULACION     
    Hora_Actual=datestr(now);
    Hora_Actual1=Hora_Actual(13:14);
    Hora_Actual2=Hora_Actual(16:17);
    Hora_Actual3=Hora_Actual(19:20);
%% NOMBRE DEL ARCHIVO A GUARDAR    
    filename=strcat(sistema, '_', Particle_name, '_',date, '_',Metaheuristica_Usada,'_',SI_continge,...
        '_',PERDAS','_','Hora','_',Hora_Actual1,'_', Hora_Actual2 ,'_',Hora_Actual3,'_','N_Etap','_',E_dinamicas,'FIN');  % Escribe un archivo .mat 

    %% DATOS ESTAD?STICOS DE LAS N PRUEBAS REALIZADAS
    SucRunStats = RunStats(RunStats(:,1)==1,:);
    ITER_PROMEDIO = mean(SucRunStats(:,2));
    ITER_DESV_St=std(SucRunStats(:,2));
    EVA_PROMEDIO_FO= mean(SucRunStats(:,3));
    DESV_St__EVA_FO= std(SucRunStats(:,3));

    %% GUARDAS RESULTADOS EN WIZARD(MATLAB)
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
                'DESV_St__EVA_FO','Tiempo_total_de_simlacion') % Graba el archivo .mat con el resultado de algunas variables

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
% Tiempo_total_de_simlacion: tiempo en el que se tarda en encontrar el resultado final 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %   for N_exp=1:nexp
 %           [xopt]=Grafica_sistema(xopt(:,N_exp),...
 %          sistema,xopt,N_exp,incluir_contingencias,Metaheuristica_Usada,incluir_perdidas);   
  %          Compensacion_actual_verificada(:,N_exp)=Comp_shunt_en_cada_barra(:,N_exp);    
 %       Gen_Pot_Activa_todas_etapas(N_exp)=sum(Gen_Pot_Activa_etapa(N_exp:N_exp)); % Calcula la generacion total de todos los estados    
 %       Costo_Gen_Pot_Activa_todas_etapas(N_exp)=sum(Costo_Gen_Pot_Activa_etapa(N_exp:N_exp));
 %   end
end

% --- Executes just before TNEP_Dinam is made visible.
function TNEP_Dinam_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TNEP_Dinam (see VARARGIN)

% Choose default command line output for TNEP_Dinam
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TNEP_Dinam wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = TNEP_Dinam_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit1_Callback(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
if str2double(get(hObject,'String'))<=0 
errordlg('Dinamic states must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Dinamic states must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('Dinamic states be a integer','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
 
   % --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(~, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit3_Callback(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Size population must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Size population must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('Size population be a integer','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, ~, ~)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit4_Callback(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
if str2double(get(hObject,'String'))<0 
  errordlg('Refer. value must be>0 (TNEP is always positive)','Error'); 
end
if str2double(get(hObject,'String'))<0 
errordlg('Refer. value must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Refer. value must be a number','Error')
end

end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, ~, ~)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit5_Callback(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Max. Itertations must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Max. Itertations must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('Max. Itertations be a integer','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, ~, ~)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit6_Callback(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Test number must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Test number must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('Test number be a integer','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, ~, ~)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(~, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(~, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit7_Callback(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Losses must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Losses must be a number','Error')
end


end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, ~, ~)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit8_Callback(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Lambda must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Lambda must be a number','Error')
end

end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, ~, ~)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit9_Callback(hObject, ~, ~)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Max. Lines must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Max. Lines must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('Max. Lines be a integer','Error');
end

end

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, ~, ~)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Increase demand must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Increase demand must be a number','Error')
end
end

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, ~, ~)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(~, ~, ~)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(~, ~, ~)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
end

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit11_Callback(hObject, ~, ~)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
if str2double(get(hObject,'String'))<0 
errordlg('Interest rate must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('Interest rate must be a number','Error')
end
end

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, ~, ~)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit12_Callback(hObject, ~, ~)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
if str2double(get(hObject,'String'))<0 
errordlg('# of iteration (repetitions) must be a number>0','Error');  
elseif isnan(str2double(get(hObject,'String'))) 
errordlg('# of iteration (repetitions) must be a number','Error')
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
  errordlg('# of iteration (repetitions) be a integer','Error');
end

end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, ~, ~)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit13_Callback(hObject, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
if str2double(get(hObject,'String'))<0 
  errordlg('Max. Eval Funtions must be>0','Error');
elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
errordlg('Max. Eval Funtions must be a integer','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%function edit15_Callback(hObject, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
%if str2double(get(hObject,'String'))<1 
%  errordlg('Initial population must be>0','Error');
%elseif  (not(mod((str2double(get(hObject,'String'))),1)))==0
%errordlg('Initial population must be a integer','Error');
%end
%end

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, ~, ~)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6
end

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7
end

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
end

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
