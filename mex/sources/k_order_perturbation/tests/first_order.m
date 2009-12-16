function [gy]=first_order(M_, dr, jacobia)
% Emulation of Dynare++ c++ first_order.cpp for testing pruposes

% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% fd = jacobia_
% reorder jacobia_
[fd]=k_reOrderedJacobia(M_, jacobia)

%ypart=dr;
ypart.ny=M_.endo_nbr;
ypart.nyss=dr.nfwrd+dr.nboth;
ypart.nys=dr.npred;
ypart.npred=dr.npred-dr.nboth;
ypart.nboth=dr.nboth;
ypart.nforw=dr.nfwrd;
ypart.nstat =dr.nstatic
nu=M_.exo_nbr

off= 1; % = 0 in C
fyplus = fd(:,off:off+ypart.nyss-1);
off= off+ypart.nyss;
fyszero= fd(:,off:off+ypart.nstat-1);
off= off+ypart.nstat;
fypzero= fd(:,off:off+ypart.npred-1);
off= off+ypart.npred;
fybzero= fd(:,off:off+ypart.nboth-1);
off= off+ypart.nboth;
fyfzero= fd(:,off:off+ypart.nforw-1);
off= off+ypart.nforw;
fymins= fd(:,off:off+ypart.nys-1);
off= off+ypart.nys;
fuzero= fd(:,off:off+nu-1);
off=off+ nu;

n= ypart.ny+ypart.nboth;
%TwoDMatrix 
matD=zeros(n,n);
%	matD.place(fypzero,0,0);
matD(1:n-ypart.nboth,1:ypart.npred)= fypzero;
%	matD.place(fybzero,0,ypart.npred);
matD(1:n-ypart.nboth,ypart.npred+1:ypart.npred+ypart.nboth)=fybzero;
%	matD.place(fyplus,0,ypart.nys()+ypart.nstat);
matD(1:n-ypart.nboth,ypart.nys+ypart.nstat+1:ypart.nys+ypart.nstat+ypart.nyss)=fyplus;
for i=1:ypart.nboth
    matD(ypart.ny()+i,ypart.npred+i)= 1.0;
end

matE=[fymins, fyszero, zeros(n-ypart.nboth,ypart.nboth), fyfzero; zeros(ypart.nboth,n)];
%	matE.place(fymins;
%	matE.place(fyszero,0,ypart.nys());
%	matE.place(fyfzero,0,ypart.nys()+ypart.nstat+ypart.nboth);

for i= 1:ypart.nboth
    matE(ypart.ny()+i,ypart.nys()+ypart.nstat+i)= -1.0;
end
matE=-matE; %matE.mult(-1.0);

%    vsl=zeros(n,n);
%	vsr=zeros(n,n);
%	lwork= 100*n+16;
%	work=zeros(1,lwork);
%	bwork=zeros(1,n);
%int info;

%    	LAPACK_dgges("N","V","S",order_eigs,&n,matE.getData().base(),&n,
%		matD.getData().base(),&n,&sdim,alphar.base(),alphai.base(),
%		beta.base(),vsl.getData().base(),&n,vsr.getData().base(),&n,
%		work.base(),&lwork,&(bwork[0]),&info);

[matE1,matD1,vsr,sdim,dr.eigval,info] = mjdgges(matE,matD,1);

bk_cond= (sdim==ypart.nys);

%  	ConstGeneralMatrix z11(vsr,0,0,ypart.nys(),ypart.nys());
z11=vsr(1:ypart.nys,1:ypart.nys);
%	ConstGeneralMatrix z12(vsr,0,ypart.nys(),ypart.nys(),n-ypart.nys());
z12=vsr(1:ypart.nys(),ypart.nys+1:end);%, n-ypart.nys);
                                       %	ConstGeneralMatrix z21(vsr,ypart.nys(),0,n-ypart.nys(),ypart.nys());
z21=vsr(ypart.nys+1:end,1:ypart.nys);
%	ConstGeneralMatrix z22(vsr,ypart.nys(),ypart.nys(),n-ypart.nys(),n-ypart.nys());
z22=vsr(ypart.nys+1:end,ypart.nys+1:end);

% 	GeneralMatrix sfder(z12,"transpose");
sfder=z12';%,"transpose");
           %	z22.multInvLeftTrans(sfder);
sfder=z22'\sfder;
sfder=-sfder;% .mult(-1);

%s11(matE,0,0,ypart.nys(),ypart.nys());
s11=matE1(1:ypart.nys,1:ypart.nys);
%	 t11=(matD1,0,0,ypart.nys(),ypart.nys());
t11=matD1(1:ypart.nys,1:ypart.nys);
dumm=(s11');%,"transpose");
            %z11.multInvLeftTrans(dumm);
dumm=z11'\dumm;
preder=(dumm');%,"transpose");
               %t11.multInvLeft(preder);
preder=t11\preder;
%preder.multLeft(z11);
preder= z11*preder;

%	gy.place(preder,ypart.nstat,0);
%	gy=(zeros(ypart.nstat,size(preder,2)) ;preder);
%	 sder(sfder,0,0,ypart.nstat,ypart.nys());
sder=sfder(1:ypart.nstat,1:ypart.nys);
%	gy.place(sder,0,0);
%	gy(1:ypart.nstat, 1:ypart.nys)=sder;
%    gy=[sder;preder];
%	 fder(sfder,ypart.nstat+ypart.nboth,0,ypart.nforw,ypart.nys());
fder=sfder(ypart.nstat+ypart.nboth+1:ypart.nstat+ypart.nboth+ypart.nforw,1:ypart.nys);
%	gy.place(fder,ypart.nstat+ypart.nys(),0);
%	gy(ypart.nstat+ypart.nys,1)=fder;
gy=[sder;preder;fder];
