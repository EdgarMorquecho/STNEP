function [poblacion, fevalcount_creacion]=creacioninicial(fhd,...
sistema,Comp_Indu_Capa_etap,Contin_Segun_Dimensi_Final,...
incluir_contingencias,incluir_perdidas,factor_perdidas,...
compensacion_anual_optim,landa,Comp_Indu_Capa_inicial,...      
SS,Xmin,Xmax,dim,opcion)
%% variables necesarias para la solucion inicial estatica            
poblacion=zeros(dim,SS);
fevalcount_creacion=0;
    switch opcion
%% 1. Método Tradicional (Aleatorio)
    case {'Aleatoria'}        
%% -Estático--
            for i=1:SS
                for j=1:dim
                    poblacion(j,i)=round(Xmin(j,1)+((Xmax(j,1))-Xmin(j,1))*(rand));
                end
            end            
        
    case {'Seudoaleatoria'}
%% -Estático--
            for i=1:SS
                for j=1:dim
                    poblacion(j,i)=round(Xmin(j,1)+((Xmax(j,1))-Xmin(j,1))*(rand));
                end
            end 
               for i=1:SS 
                   indices=randperm(dim,round(0.65*dim));
                   poblacion(indices,i)=Xmin(indices,1);
               end
                              
    
    end
end
