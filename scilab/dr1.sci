function [dr]=dr1(iorder,dr,correct)
// Copyright (C) 2001 Michel Juillard
// 
 
xlen = xkmax_+xkmin_+1;
klen = ykmin_+ykmax_+1;
iyv = iy_';
iyv = iyv(:);
iyr0 = find(iyv)';
it_ = ykmin_+1;
 
 
if exo_nbr==0 then
  exe_ = [];
end
 
if ~(iy_(ykmin_+1,:)>0) then
  error('Error in model specification: some variables don""t appear as current');
end
 
if ykmax_==0 then
  error('No forward variable: no point in using dr1');
end
 
if xlen>1 then
  error('SS: stochastic exogenous variables must appear only at the'+' current period. Use additional endogenous variables');
end
 
tempex = ex_;
 
it_ = ykmin_+1;
 
z = dr('ys')*ones(1,klen);
z = z(iyr0);
 
jacobia_ = real(jacob_a('ff1_',[z;exe_]));
 
ex_ = tempex;
tempex = [];
 
nz = size(z,1);
pred_var = find(or(iy_(1:ykmin_,:),1));
 
fwrd_var = find(or(iy_(ykmin_+2:ykmin_+ykmax_+1,:),1));
both_var = intersect(pred_var,fwrd_var);
if size(both_var,'*')
  pred_var = setdiff(pred_var,both_var);
  fwrd_var = setdiff(fwrd_var,both_var);
end
nboth = length(both_var);
npred = length(pred_var);
nfwrd = length(fwrd_var);
stat_var = setdiff(1:endo_nbr,union(union(pred_var,both_var),fwrd_var));
nstatic = length(stat_var);
order_var = [stat_var,pred_var,both_var,fwrd_var];
k1 = iy_(find((1:klen)~=(ykmin_+1)),:);
b = jacobia_(:,iy_(ykmin_+1,order_var));
a = b\jacobia_(:,nonzeros(k1'));
if exo_nbr then
  fu = b\jacobia_(:,nz+1:size(jacobia_,2));
end
 
// building kmask for z state vector in t+1
if ykmin_>0 then
  if ykmax_>0 then
    %v2 = iy_(ykmin_+2:size(iy_,1),order_var)
    kmask = cumsum(%v2($:-1:1,:),1);
    %v2 = cumsum(iy_(1:ykmin_,order_var),1)
    kmask = [kmask;%v2($:-1:1,:)];
  end
else
  %v1 = iy_(ykmin_+1:klen,order_var)
  kmask = cumsum(%v1($:-1:1,:),1);
end
 
kmask = kmask';
kmask = kmask(:);
i_kmask = find(kmask)';
// index of nonzero entries in kmask
nd = size(i_kmask,1);
// size of the state vector
%v = 1:nd
kmask(i_kmask,1) = %v(:);
 
// auxiliary equations
// elements that are both in z(t+1) and z(t)
k1 = find(kmask(1:length(kmask)-endo_nbr)&kmask(endo_nbr+1:length(kmask)))';
kad = [];
kae = [];
if ~(k1==[]) then
  kad = kmask(k1+endo_nbr);
  kae = kmask(k1);
end
 
// composition of state vector
// col 1: variable;           col 2: lead/lag in z(t+1); 
// col 3: A cols for t+1 (D); col 4: A cols for t (E)
 
kstate = [ones(klen-1,1).*.[1:endo_nbr]' ((klen:-1:2)').*.ones(endo_nbr,1) zeros((klen-1)*endo_nbr,2)];
%v = iy_(:,order_var)
kiy = %v($:-1:1,:)';
kiy = kiy(:);
kstate(1:ykmax_*endo_nbr,3) = kiy(1:ykmax_*endo_nbr)-endo_nbr;
kstate(ykmax_*endo_nbr+1:size(kstate,1),4) = kiy((ykmax_+1)*endo_nbr+1:length(kiy));
// put in E only the current variables that are not already in D
kstate = kstate(i_kmask,:);
dr('kstate') = kstate;

sdyn = endo_nbr-nstatic;
 
// buildind D and E
d = zeros(nd,nd);
e = d;
 
k = find((kstate(:,2)>=(ykmin_+2)) & kstate(:,3))';
d(1:sdyn,k) = a(nstatic+1:size(a,1),kstate(k,3));
k1 = find(kstate(:,2)==(ykmin_+2))';
 
a1 = eye(sdyn,sdyn);
e(1:sdyn,1:size(k1,1)) = -a1(:,kstate(k1,1)-nstatic);
k = find((kstate(:,2)<=(ykmin_+1))&kstate(:,4))';
e(1:sdyn,k) = -a(nstatic+1:size(a,1),kstate(k,4));
k2 = find(kstate(:,2)==(ykmin_+1))';
[k3,k4] = setdiff(kstate(k2,1),kstate(k1,1));

if size(k3,'*') > 0 then
  d(1:sdyn,k2(k4)) = a1(:,k3-nstatic);
end

if ~(kad==[]) then
  for j = 1:size(kad,1)
    d(sdyn+j,kad(j)) = 1;
    e(sdyn+j,kae(j)) = 1;
  end
end

if grep(getversion(),'2.6') then
  [ss,tt,w,sdim] = gschur(e,d,'d');
else
  [ss,tt,w,sdim] = schur(e,d,'d');
end  

nba = nd-sdim;
nyf = sum(kstate(:,2)>(ykmin_+1));
if nba~=nyf then
  dyn_disp('LSS warning: Blanchard-Kahn conditions are not satisfied. Run CHECK to learn more!');
  halt();
   
end
 
np = nd-nyf;
n2 = np+1;
n3 = nyf;
n4 = n3+1;
// derivatives with respect to dynamic state variables
// forward variables
gx = -w(1:n3,n2:nd)'\w(n4:nd,n2:nd)';
// predetermined variables
hx = w(1:n3,1:np)'*gx+w(n4:nd,1:np)';
hx = tt(1:np,1:np)*hx\(ss(1:np,1:np)*hx);
 
 
k1 = find(kstate(n4:nd,2) == ykmin_+1);
k2 = find(kstate(1:n3,2) == ykmin_+2);
dr('ghx') = [hx(k1,:);gx(k2(nboth+1:$),:)]
 
// lead variables actually present in the model
j3 = nonzeros(kstate(:,3));
j4  = find(kstate(:,3));
// derivatives with respect to exogenous variables
if exo_nbr then
  a1 = eye(endo_nbr,endo_nbr);
  aa1 = [];
  if nstatic>0 then
    aa1 = a1(:,1:nstatic);
  end
  dr('ghu') = -[aa1,a(:,j3)*gx(j4,1:npred+nboth)+a1(:,nstatic+1:nstatic+npred+nboth),a1(:,nstatic+npred+nboth+1:$)]\fu
end
 
// static variables
if nstatic>0 then
    temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
    j5 = find(kstate(n4:nd,4));
    temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
    dr.ghx = [temp; dr.ghx];
    temp = [];
end
 
dr('order_var') = order_var
dr('nstatic') = nstatic
dr('npred') = npred+nboth
 
if iorder==1 then
  return
end
 
// Second order
tempex = ex_;
 
hessian = real(jacob2('ff1_',[z;exe_]));

ex_ = tempex;
clear('tempex');
 
 
np1 = sum(sum(iy_(1:ykmin_,:)>0,2));
zx = [eye(npred+nboth,npred+nboth);dr('ghx');gx*hx;zeros(exo_nbr,np)];
 k1 = nonzeros(iy_(:,order_var)');
if exo_nbr then
  k1 = [k1;(length(k1)+1:length(k1)+exo_nbr)'];
end
kk = matrix(1:length(k1)^2,length(k1),length(k1));
// reorder hessian
kk = kk(k1,k1);
hessian = hessian(:,kk);
rhs = -hessian*(zx.*.zx);
rhs = rhs(:);
 
//lhs
hxt = hx';
gx = dr('ghx');
gx = gx(size(gx,1)-nyf+1:size(gx,1),:); 
Inp = speye(np*np,np*np);
 
f1 = sparse(jacobia_(:,nonzeros(iy_(ykmin_+2,:))));
// f_y'g_xx hx kron hx + f_y g_xx + f_y'g_x h_xx + f_x' h_xx
// 
j1 = sparse(jacobia_(:,iy_(ykmin_+1,order_var)));
lhs = Inp.*.j1;
k = matrix(1:endo_nbr*np*np,endo_nbr,np*np);
nr = size(k,1);
lhs(:,k(nr-nyf+1:nr,:)) = lhs(:,k(nr-nyf+1:nr,:))+hxt.*.hxt.*.f1;
lhs(:,k(nstatic+1:nstatic+np,:)) = lhs(:,k(nstatic+1:nstatic+np,:))+Inp.*.(f1*gx);
dr('ghxx') = lhs\rhs
dr('ghxx') = matrix(dr('ghxx'),endo_nbr,np*np)
 
//ghxu
//rhs
hu = dr('ghu')(nstatic+1:nstatic+np,:);
zu = [zeros(np1,exo_nbr);dr('ghu');gx*hu;eye(exo_nbr,exo_nbr)];
nr = size(dr('ghxx'),1);
rhs = -hessian*(zx.*.zu)-f1*dr('ghxx')(nr-nyf+1:nr,:)*(hx.*.hu);
 
//lhs
lhs = j1;
lhs(:,nstatic+1:nstatic+np) = lhs(:,nstatic+1:nstatic+np)+f1*gx;
dr('ghxu') = lhs\rhs
 
//ghuu
//rhs
rhs = -hessian*(zu.*.zu)-f1*dr('ghxx')(nr-nyf+1:nr,:)*(hu.*.hu);
 
//lhs
dr('ghuu') = lhs\rhs
 
 
// dr.ghs2
// derivatives of F with respect to forward variables
k1 = matrix(1:size(hessian,2),size(jacobia_,2),size(jacobia_,2));
k1 = k1(npred+nboth+endo_nbr+1:npred+endo_nbr+2*nboth+nfwrd,npred+nboth+endo_nbr+1:npred+endo_nbr+2*nboth+nfwrd);
 
nr = size(dr('ghuu'),1);
guu = dr('ghuu')(nr-nyf+1:nr,:);
nr = size(dr('ghu'),1);
gu = dr('ghu')(nr-nyf+1:nr,:);
rhs = -(f1*guu+hessian(:,k1)*(gu.*.gu))*Sigma_e_(:);
lhs = j1;
lhs(:,nstatic+1:nstatic+np) = lhs(:,nstatic+1:nstatic+np)+f1*gx;
nc = size(lhs,2);
lhs(:,nc-nyf+1:nc) = lhs(:,nc-nyf+1:nc)+f1;
dr('ghs2') = lhs\rhs;
 
// 01/08/2001 MJ put return on iorder == 1 after defining dr.kstate and dr.kdyn
// 01/17/2001 MJ added dr.delta_s: correction factor for order = 2
// 01/21/2001 FC correction of delta_s for more than 1 shock
// 01/23/2001 MJ suppression of redundant sum() in delta_s formula
// 02/22/2001 MJ stderr_ replaced by Sigma_e_
// 08/24/2001 MJ changed the order of the variables, separates static
//               variables and handles only one instance of variables both
//               in lead and lag
// 08/24/2001 MJ added sigma to state variables as in Schmitt-Grohe and
//               Uribe (2001)
// 10/30/2002 MJ corrected bugs with lags on several periods 
 
 
 
 
 
 
 
 
 
