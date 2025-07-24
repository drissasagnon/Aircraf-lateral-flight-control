function [MOD,COR] = clickmap(MOD_IN,COR_IN,str,bid);
%
% [MOD_OUT,COR_OUT] = CLICKMAP(MOD_IN,COR_IN)
%
% � partir du lieu des racines et de mani�re interactive, 
% CLICKMAP permet de s�lectionner et de supprimer
% des modes du mod�le MOD_IN et/ou du correcteur COR_IN (dans leurs 
% bases modales r�elles respectives) et de redessiner le lieu des 
% racines correspondant pour valider ou non la s�lection.
%
% MOD_IN et COR_IN sont les matrices syst�mes ou les objets LTI repr�sentants
% respectivement le mod�le et le correcteur.
%
% MOD_OUT et COR_OUT sont les matrices syst�mes repr�sentants le mod�le 
% et le correcteur r�duits dans leurs bases diagonales.
%
% ATTENTION!: - les bases modales r�elles sont obtenues par la fonction
%               CANON. Il peut y avoir des probl�mes en cas de valeurs
%               propres multiples (bloc de Jordan).
%
%             - Le rebouclage du correcteur sur le mod�le est suppos�
%               positif.
%
% [MOD_OUT,COR_OUT] = CLICKMAP(MOD_IN,COR_IN,'s') permet de travailler
% avec des r�alisations de Schur pour le mod�le et le correcteur (on
% tronque dans les formes de Schur r�elle: pr�f�rable en cas de valeurs
% propres multiples).
%
%
% Voir aussi: RLOCUSP, CLICKPOLE

%  D. Alazard 02/99
%  Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.

global niveau
if isempty(niveau), niveau=0;end;

global KK
global mygcf

if nargin==2,str='m';end;
if nargin<=3,bid=0;end;

[a,b,c,d,flag,T]=test_natss(MOD_IN);

[ak,bk,ck,dk,flagk,Tk]=test_natss(COR_IN);

valpk=eig(ak);
ck=-ck;dk=-dk;

if bid~=1,
        [ao,bo,co,do]=series(a,b,c,d,ak,bk,ck,dk);
	if bid==0,
	   valpbo=eig(ao);
      ax=max(max(abs(real(valpbo))),1);
      ay=max(max(abs(imag(valpbo))),1);
      ax=max(ax,ay);ax=[-1.2*ax 1.2*ax];
      ay=ax;
	else
	   ax=get(gca,'Xlim');
	   ay=get(gca,'Ylim');
	end
	figure('Name',['r�duction niveau: ',num2str(niveau)]);
	portrait
	axis([ax ay]);
	hold on
   rlocusp(ao,bo,co,do);
   plot(real(valpk),imag(valpk),'gx','LineWidth',2);
	if T==0, sgrid, else, zgrid;end;
	plot([0 0],ay);
	plot(ax,[0 0]);
end	
	
if niveau==0, mygcf=gcf;end;
KK=1;
while KK==1,
  disp(' ');
  disp(' ');
  disp(['Etape niveau: ',num2str(niveau)]);
  disp(' ');
  disp(['Ordre du syst�me : ',num2str(size(a,1))]);
  disp(['Ordre du correcteur : ',num2str(size(ak,1))]);
  figure(mygcf+niveau)
  zoom on
  if niveau==0,
      KK=menu('Reduction',...
                      'd�marrer �tape 1',...
                      'stop');
  else
      KK=menu('Reduction',...
                     ['d�marrer etape ',num2str(niveau+1)],...
                     ['annuler etape ',num2str(niveau)],'stop');
  end
  if KK==1,
      figure(mygcf+niveau)
      MOD_1=MOD_IN;
      COR_1=COR_IN;
      [MOD,COR,indic1,indic2] = clickpole(MOD_1,COR_1,str);
      if indic1~=0|indic2~=0,
         niveau=niveau+1;
         if (nargout==1)&(indic2~=0),
            disp('============================================================');
            disp(' Attention: le correcteur a �t� r�duit mais il n''est pas en');
            disp(' argument de sortie de la fonction');
         end
         MOD_1=MOD;
         COR_1=COR;
         if nargout==2, 
             [MOD,COR]=clickmap(MOD_1,COR_1,str,2);
         else
             MOD=clickmap(MOD_1,COR_1,str,2);
         end
         if  KK==2,
            delete(mygcf+niveau)
            niveau=niveau-1;
            MOD=MOD_IN;
            if nargout==2, COR=COR_IN;end;
            KK=1;
        end
      end
  else
      MOD=MOD_IN;
      if nargout==2, COR=COR_IN; end;
      if KK==3, clear global niveau mygcf,end;
  end
end


function [as,bs,cs,ds,aks,bks,cks,dks,indic,indick]=clickpole(a,b,str)

[ak,bk,ck,dk,flagk,Tk] = test_natss(b);
[a,b,c,d,flag,T] = test_natss(a);
if str=='s',
   [U,ak]=schur(ak);
   bk=U'*bk,ck=ck*U;
   [U,a]=schur(a);
   b=U'*b,c=c*U;
else
   [ak,bk,ck,dk] = canon(ak,bk,ck,dk,'modal');
   [a,b,c,d] = canon(a,b,c,d,'modal');
end;	

indic=0;
indick=0;

na = size(a,1); nak = size(ak,1);
matrice = [a zeros(na,nak);zeros(nak,na) ak];

index = polefind2(matrice);

liste_a = find(index <= na);
liste_k = find(index >= na+1);
index_a = index(liste_a);
index_k = index(liste_k);
index_k = index_k - na;

as = a; bs = b; cs = c; ds = d; aks = ak; bks = bk; cks = ck; dks = dk;

if ~isempty(index_a), 

[as,bs,cs,ds]=modred(a,b,c,d,index_a);
disp(' ');
disp(['Vous avez r�duit ',num2str(length(index_a)),' p�le(s) dans le syst�me.']);
indic=length(index_a);
end;

if ~isempty(index_k),

[aks,bks,cks,dks]=modred(ak,bk,ck,dk,index_k);
disp(' ');
disp(['Vous avez r�duit ',num2str(length(index_k)),' p�le(s) dans le correcteur.']);
indick=length(index_k);
end;

if nargout <= 4,
	if flag == 2,
	                as=pck(as,bs,cs,ds);
	else, 
	  		as=ss(as,bs,cs,ds,T);
	end;		
	if flagk == 2,   
			bs=pck(aks,bks,cks,dks);
	else,  		
			bs=ss(aks,bks,cks,dks,Tk);
	end;		
	cs=indic;ds=indick;	
end;

function index=polefind2(a,nb)
% INDEX = POLEFIND(A,NB) permet de selectionner a l'aide de la souris 
% NB modes de A sur un lieu des racines ou une carte des poles de A
% pr�alablement trac�.
% INDEX retourne les indices de A correspondant aux modes selectionnes
%
% IMPORTANT: A doit etre sous forme Bloc Diagonale Reelle
% 
% Un pole reel correspond a une ligne de A
% Un pole complexe s�lectionn� , et obligatoirement son conjugue,
% correspondent ensemble a deux lignes successives de A { bloc matriciel 2*2 }
%
% Si NB est specifie, la fonction attend le cliquage de NB poles pour remplir
% le vecteur INDEX. (Quand un pole complexe est selectionne, son conjugue
% l'est aussi, donc LENGTH(INDEX) peut etre superieure a NB )
% 
% Si NB n'est pas specifie, le remplissage de INDEX s'effectue au coup par
% coup, par cliquage successif de pole. Les boutons de menu apparaissent alors
% sur la barre d'outil de la figure courante.
%
% SEE ALSO  EIG_BDIAG, RBDIAG
%

% S. Delonnoy 01/99
% Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.
%

index = [];
zoom on
hold on
if nargin ==2,

V=eig_bdiag(a);
taille=size(V,1);
reV=real(V);
imV=imag(V);

if nb == 1, fprintf(['cliquez le mode que vous desirez\n']);
else fprintf(['cliquez les ',num2str(nb),' modes que vous desirez\n']);
end;

liste=ginput(nb);

i=sqrt(-1);
u=1;
for k=1:nb,
			test1 = V-(liste(k,1)+i*liste(k,2))*ones(taille,1);
			abs_test1 = abs(test1);
			[mini1,ind] = min(abs_test1);
			message = redondance(ind,index,k);
		if message ~=1,
			index(u) = ind;
			u = u+1;
		if abs(imV(ind)) > 1e-6*abs(reV(ind)),
			test2 = V-(liste(k,1)-i*liste(k,2))*ones(taille,1);
			abs_test2 = abs(test2);
			[mini2,ind2] = min(abs_test2);
			ind1 = min(ind,ind2);
			ind2 = max(ind,ind2);

			index(u-1) = ind1;
			index(u) = ind2;
			u = u+1;
		end;
		end;

end;
index = index';

elseif nargin ==1,

KK = 1;
%disp(' ');
%disp('s�lectionner les modes que vous desirez en cliquant d''abord sur "s�lection",')
%disp('corriger �ventuellement la s�lection en cliquant sur "correction",')
%disp('arr�ter en cliquant sur "STOP" ');

index=[];
handler=[];
%hmenu = mon_menu(['Faites votre choix';'   (ZOOM ACTIF)   '],'s�lection de modes','d�selection','stop');
while KK ~= 3
   %        KK=mon_menu(hmenu);
   KK=menu(['Faites votre choix';'   (ZOOM ACTIF)   '],'s�lection de modes','d�selection','stop');
	if KK == 1, ind=nouveau_point(index,a);
	     message = redondance(ind(1),index);
	     if ~message,
	            handi=plot(real(eig(a(ind,ind))),imag(eig(a(ind,ind))),'c*','LineWidth',2);
	            index=[index;[ind]];
	            handler=[handler;handi];
	            if length(ind)==2,
	               handler=[handler;handi];
	            end;
	     end
	elseif KK == 2, ind=nouveau_point(index,a);
	     if ~isempty(ind),
	            indice=find(index==ind(1));
	            delete(handler(indice));
	            indice=find(index~=ind(1));
	            if length(ind)==2,
	               indice=find((index~=ind(1))&(index~=ind(2)));
	            end;
	            if ~isempty(handler) handler=handler(indice);end;
	            if ~isempty(index) index=index(indice);end;
	     end
	end;
end;
% On efface les marqueurs
delete(handler)
%delete(hmenu)
% On retourne � la pleine �chelle
zoom out
else, 
error('les arguments d''entree sont une matrice et eventuellement un nombre')
return
end;

% =============================================
% fonction interne
% =============================================

function indo=nouveau_point(index,a)
% rajoute un mode a la selection

indo=[];
V=eig_bdiag(a);
taille=size(V,1);
reV=real(V);
imV=imag(V);

	liste = ginput(1);
	
	i=sqrt(-1);
			test1 = V-(liste(1)+i*liste(2))*ones(taille,1);
			abs_test1 = abs(test1);
			[mini1,ind] = min(abs_test1);
			indo=ind;
		if abs(imV(ind)) > 1e-6*abs(reV(ind))
			test2 = V-(liste(1)-i*liste(2))*ones(taille,1);
			abs_test2 = abs(test2);
			[mini2,ind2] = min(abs_test2);
			ind1 = min(ind,ind2);
			ind2 = max(ind,ind2);

			indo=[ind1;ind2];
			
		end;

% end of function nouveau_point


% ==========================================
% fonction interne
% ==========================================
function message=redondance(ind,index,k)
% verifie que le mode n'a pas deja ete selectionne
%
message = 0;
if ~isempty(index),
xx = find(index==ind);
if ~isempty(xx), 
		if nargin == 2,
                disp('============================================ ');
		disp('ce mode a deja ete selectionne');
  		elseif nargin == 3,
                disp('============================================ ');
fprintf(['le ',num2str(k),'-ieme mode choisi a deja ete selectionne\n']);
  		end;
  		message = 1;
end; 
end;
% end of function redondance

function [vec]=eig_bdiag(a)
%
% [VEC]=EIG_BDIAG(A) trouve les valeurs propres de la matrice 
% bloc Diagonale Reelle A dans l'ordre de leur 'apparition' dans A.
% VEC contient donc les valeurs propres des blocs successifs constituant A
%
% ATTENTION: A obligatoirement sous forme Bloc Diagonale Reelle
%

% N. Imbert /01/01/96
% Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.

%
s=size(a,1);
j=1;
vec=zeros(s,1);

while j<s,
	if abs(a(j,j+1))<=1e-6*abs(a(j,j)) & abs(a(j+1,j))<=1e-6*abs(a(j,j)),
		vec(j)=a(j,j);
		j=j+1;
	
	else    b=a(j:j+1,j:j+1);
		vb=eig(b);
		for k=[0 1],vec(j+k)=vb(1+k);
		end;
		j=j+2;
	end;
end;
if j==s, vec(s)=a(s,s);
end;

		