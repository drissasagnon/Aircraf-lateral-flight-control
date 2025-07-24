function [sysout]=red_fast(sysin,w)
% [SYS_OUT]=RED_FAST(SYS_IN,W)
%  enl�ve dans le syst�me SYS_IN (object LTI) les modes stables
%  et rapides dont la partie r�elle est inf�rieure � W
%  et conserve le gain statique.
%
%  Repr�sentation du syst�me r�duit syst�me LTI (ss) SYS_OUT.
%     

% D. Alazard 01/95
% Copyright (c) 1993-2000 ONERA/DCSD, All Rights Reserved.
% Revised 03/2011

sys=canon(sysin,'modal');

[aa,bb,cc,dd]=ssdata(sys);
n=size(aa,1);
as=[];bs=[];cs=[];ds=0*dd;
gain=dcgain(sysin);
j=1;
i=1;
while i<=n-1,
    if aa(i,i+1)==0,  % real eigenvalue
        if real(aa(i,i))>w,
            as(j,j)=aa(i,i);
            bs(j,:)=bb(i,:);
            cs(:,j)=cc(:,i);
            j=j+1;
        else,
    		gain_i=-cc(:,i)/aa(i,i)*bb(i,:);
    		ds=ds+gain_i;
        end;
        i=i+1;
    else % complex eigenvalue
        if real(aa(i,i))>w,
            as(j:j+1,j:j+1)=aa(i:i+1,i:i+1);
            bs(j:j+1,:)=bb(i:i+1,:);
            cs(:,j:j+1)=cc(:,i:i+1);
            j=j+2;
        else,
    		gain_i=-cc(:,i:i+1)*inv(aa(i:i+1,i:i+1))*bb(i:i+1,:);
    		ds=ds+gain_i;
        end;
        i=i+2;      
    end
end;
if i==n,  % real eigenvalue
    if real(aa(i,i))>w,
        as(j,j)=aa(i,i);
        bs(j,:)=bb(i,:);
        cs(:,j)=cc(:,i);
        j=j+1;
    else,
        gain_i=-cc(:,i)/aa(i,i)*bb(i,:);
        ds=ds+gain_i;
    end;
end
ds=ds+dd;
sysout=ss(as,bs,cs,ds);

