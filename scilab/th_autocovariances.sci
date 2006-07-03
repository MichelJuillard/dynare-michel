function [Gamma_y]=th_autocovariances(dr,ivar)
Gamma_y=[];
// Copyright (C) 2001 Michel Juillard
// 
// computes the theoretical auto-covariances, Gamma_y, for an AR(p) process 
// with coefficients dr.ghx and dr.ghu and shock variances Sigma_e_
// for a subset of variables ivar (indices in lgy_)
// 
global('lgy_','endo_nbr','exo_nbr','Sigma_e_','ykmin_');
 
nvar = size(ivar,1);
 
ghx = dr('ghx');
ghu = dr('ghu');
npred = dr('npred');
nstatic = dr('nstatic');
kstate = dr('kstate');
order = dr('order_var');
iv(1,order) = 1:endo_nbr;
nx = size(ghx,2);
 
ikx = nstatic+1:nstatic+npred;
 
A = zeros(nx,nx);
A(1:npred,:) = ghx(ikx,:);
offset_r = npred;
offset_c = 0;
i0 = find(kstate(:,2)==(ykmin_+1))';
n0 = size(i0,1);
for i = ykmin_:-1:2
  i1 = find(kstate(:,2)==i)';
  n1 = size(i1,1);
  j = zeros(n1,1);
  for j1 = 1:n1
    %v2 = find(kstate(i0,1)==kstate(i1(j1),1))'
    j(j1,1) = %v2(:);
  end
  A(offset_r+1:offset_r+n1,offset_c+j) = eye(n1,n1);
  offset_r = offset_r+n1;
  offset_c = offset_c+n0;
  i0 = i1;
  n0 = n1;
end
kron_A = A.*.A;
ghu1 = [ghu(ikx,:);zeros(nx-npred,exo_nbr)];
kron_ghu = ghu1.*.ghu1;
 
// vx such that vec(sigma_y) = vx * vec(Sigma_e_) (predetermined vars) 
vx = (eye(nx*nx,nx*nx)-kron_A)\kron_ghu;
kron_ghx = [];
kron_ghu = [];
 
// order of variables with preset variances in ghx and ghu
 
//! mtlb_e(iv,ivar) may be replaced by iv(ivar)
//!    iv(ivar) if iv is a vector,
//!    iv(ivar).' if iv is a matrix.
iky = mtlb_e(iv,ivar);
 
kron_ghx = ghx(iky,:).*.ghx(iky,:);
kron_ghu = ghu(iky,:).*.ghu(iky,:);
Gamma_y = matrix((kron_ghx*vx+kron_ghu)*Sigma_e_(:),nvar,nvar);
 
// 10/18/2002 MJ
 
 
 
 
