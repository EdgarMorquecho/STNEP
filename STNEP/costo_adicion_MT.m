function [costo_total]=costo_adicion_MT(x_posicion, sistema,dim,Costos)

    
%% Lineas de trasnmision adicionadas
    LT_nuevas=x_posicion(1:dim,:); % recoje los costos y la 

%% Calcular el costo por lineas adicionadas
    costo_total=sum(Costos.*LT_nuevas);

