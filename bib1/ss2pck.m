function [pack]=ss2pck(sys)
% SYSO=SS2PCK(SYSI) permet de convertir un syst�me SYSI
%  repr�sent� par un objet LTI en une matrice syst�me SYSO
%  (cr��e par PCK).
%
% Voir aussi: PCK2SS, PCK et UNPCK.


% S. Delannoy 10/98
% Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.

[ass,bss,css,dss]=ssdata(sys);
pack=pck(ass,bss,css,dss);
