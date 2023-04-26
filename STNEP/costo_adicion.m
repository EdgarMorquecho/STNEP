function [costo_total]=costo_adicion(x_posicion, sistema)

%% Esta funcion implementa la evaluacion de cada topologia de transmision
% y está basado en matpower





%% Adecuación de los parámetros para el cálculo de las restricciones
  [LineasIniciales Costos]=TestSystems(sistema);


%% Adecuación de los parámetros para el cálculo de las restricciones
numero_ramas=numel(LineasIniciales); 
%%LineasIniciales
mpc=loadcase(sistema);
LT_nuevas=x_posicion(1:numero_ramas,:)-LineasIniciales';; % recoje los costos y la 


   
  %% EVALUACION DE LA x_posicion %%
% costo_total=upso_fcn(Costos', LT_nuevas);
costo_total=sum(Costos'.*LT_nuevas);%+penalizacion;
% clear mpc
