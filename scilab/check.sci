function [result]=check()
result=[];
// Copyright (C) 2001 Michel Juillard
// 
// 
 
global ex_ it_ jacobia_

xlen = xkmax_+xkmin_+1;
[klen,ny] = size(iy_);
iyr0 = matrix(find(matrix(iy_',(ykmin_+ykmax_+1)*ny,1)),1,-1);
it_ = ykmin_+1;
tempex = ex_;
 
if valf_==0 then
  ys = ones(ykmin_+ykmax_+1,1).*.ys_;
  z = ys(iyr0);
  if exo_nbr > 0 then
    ex_ = ones(xlen,1)*exe_';
  end
  i = evstr(fname_+'_ff(z)');
else
  t = y_(:,1:ykmin_+ykmax_+1);
  z = matrix(t,size(t,1)*size(t,2),1);
  z = z(iyr0);
end
 
if ~(iy_(ykmin_+1,:)>0) then
  error('Model specification error: some variable doesn''t appear in current period.');
  
end
 
if (ykmax_==0)&(ykmin_==0) then
     error('Static model: no eigenvalue to compute.');
end

it_ = ykmin_+1;
jacob(fname_+'_ff',z);
ex_ = tempex;
clear('tempex');
 
k1 = iy_(find((((1:ykmin_+ykmax_+1)==(ykmin_+1))==%f)')',:);
temp = matrix(k1',size(k1,1)*size(k1,2),1);
k2 = temp(matrix(find(temp),1,-1));
k3 = iy_(ykmin_+1,:)';
a = jacobia_(:,k2);
b = jacobia_(:,k3);
a = b\a;
clear('b','jacobia_');
 
if ykmin_>0 then
  k2 = cumsum(iy_(1:ykmin_,:),1);
  if ykmax_>0 then
    %v2 = iy_(ykmin_+2:klen,:)
    temp = %v2($:-1:1,:);
    if size(temp,1)>1 then
      temp = cumsum(temp,1);
    end
    k2 = [k2;temp($:-1:1,:)];
  end
else
  %v1 = iy_(ykmin_+2:klen,:)
  temp = %v1($:-1:1,:);
  if size(temp,1)>1 then
    temp = cumsum(temp,1);
  end
  k2 = temp($:-1:1,:);
end
 
%temp = k2';
k3 = %temp(:);
k4 = [find(k3)',k3(find(k3))];
k3(k4(:,1)) = [1:size(k4,1)]';
kad = zeros(0,0);
 
if ny<size(k3,1) then
  k4 = k3(1:size(k3,1)-ny,:)&k3(1+ny:size(k3,1),:);
  temp = [k3(1:size(k3,1)-ny,:),k4];
  kad = temp(find(temp(:,size(temp,2)))',:);
  if nnz(kad)==size(kad,1)*size(kad,2) then
    kad = kad(:,1);
    temp = [k3(1+ny:size(k3,1),:),k4];
    kae = temp(find(temp(:,size(temp,2)))',:);
    kae = kae(:,1);
  end
end
clear('temp','k4');
 
if ykmax_>0 then
  temp = [(1:ny)',k2([ykmin_,ykmin_+1],:)'];
  k3 = temp(matrix(find(sum(k2([ykmin_,ykmin_+1],:),1)),1,-1)',:);
else
  temp = [(1:ny)',k2(ykmin_,:)'];
  k3 = temp(matrix(find(k2(ykmin_,:)'),1,-1),:);
end
 
temp = matrix(k1',size(k1,1)*size(k1,2),1);
temp = [temp,ones(klen-1,1).*.((1:ny)')];
temp = [temp,((2:klen)').*.ones(ny,1)];
k1 = temp(matrix(find(matrix(k2',size(k2,1)*size(k2,2),1)),1,-1),:);
clear('temp');
 
nd = size(k1,1);
d = zeros(nd,nd);
e = d;
 
if ykmin_>0 then
  temp = [(1:nd)',k1(:,1)];
  temp1 = (k1(:,3)<=(ykmin_+1))&(k1(:,1)>0);
  k2 = temp(find(temp1)',:);
  e(1:size(k3,1),k2(:,1)) = -a(k3(:,1),k2(:,2));
end
 
nfwrd = 0;
if ykmax_>0 then
  temp = [(1:nd)',k1(:,1)];
  temp1 = (k1(:,3)>(ykmin_+1))&(k1(:,1)>0);
  k2 = temp(find(temp1)',:);
  d(1:size(k3,1),k2(:,1)) = a(k3(:,1),k2(:,2)-ny);
  nfwrd = sum(k1(:,3) > ykmin_+1);
end
 
 
if (ykmin_>0)&(ykmax_>0) then
  a = eye(size(k3,1),size(k3,1));
  temp = (1:nd)';
  temp1 = (1:size(k3,1))';
  k2 = [temp(find(k1(:,3)==(ykmin_+2))',:),temp1(find(k3(:,3))',:)];
  e(1:size(k3,1),k2(:,1)) = -a(:,k2(:,2));
  temp1 = [(1:size(k3,1))',k3(:,3)];
  k2 = [temp(find(k1(:,3)==(ykmin_+1))',:),temp1(find(k3(:,2))',:)];
  temp = k2(:,[1,2]);
  k2 = temp(find(k2(:,3)==0)',:);
  d(1:size(k3,1),k2(:,1)) = a(:,k2(:,2));
  clear('a');
end
 
clear('temp','temp1');
 
if nnz(kad)==size(kad,1)*size(kad,2) then
  for j = 1:size(kad,1)
    d(size(k3,1)+j,kad(j)) = 1;
    e(size(k3,1)+j,kae(j)) = 1;
  end
end
 
if grep(getversion(),'2.6') then
  [ss,tt] = gschur(e,d);
else
  [ss,tt] = schur(e,d);
end

ss = diag(ss);
tt = diag(tt);

mod = ieee();
ieee(2);

lambda = ss ./ tt;

ieee(mod); 

%v = abs(lambda)
[m_lambda,%i]=sort(%v);
dyn_disp('Eigenvalues:');
dyn_disp('   Modulus    Real    Imaginary');
 
mprintf('%12.6f %12.6f\n',[m_lambda lambda(%i)]);
na = sum(abs(lambda)>1);
mprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ',na);
mprintf('for %d forward-looking variable(s)',nfwrd);
mfprintf(fh_log,'%12.6f %12.6f\n',[m_lambda lambda(%i)]);
mfprintf(fh_log,'\nThere are %d eigenvalue(s) larger than 1 in modulus ',na);
mfprintf(fh_log,'for %d forward-looking variable(s)',nfwrd);
 
// 02/06/02 MJ port to Scilab 
 
// 2/9/99 MJ: line 15, added test for absence of exogenous variable.
// 8/27/2000 MJ: change JACOB call. Added ...,1 to cumsum()
// 6/24/01 MJ: added count of abs(eigenvalues) > 1
// 02/21/03 MJ: redone output to logfile 
 
 
 
 
 
 
 
