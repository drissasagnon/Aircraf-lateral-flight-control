function [sys]=pck2ss(pack,T)
% [SYS]=PCK2SS(PACK)  permet de convertir un syst�me continu
% repr�sent� par une matrice syst�me PACK (cr��e par PCK) 
% en objet LTI d�crit par une repr�sentation d'�tat SYS.
%				                            		
% [SYS]=PCK2SS(PACK,T)  permet de convertir un syst�me �chantillonn�
% � la cadence T et repr�sent� par une matrice syst�me PACK 
% (cr��e par PCK) en objet LTI d�crit par une repr�sentation d'�tat SYS.
%				                            		
% Voir aussi: SS2PCK, PCK et UNPCK

% S. Delonnoy 01/99
% Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.

[ass,bss,css,dss]=unpck(pack);
if nargin==1,
   sys=ss(ass,bss,css,dss);
else
   sys=ss(ass,bss,css,dss,T);
end;
