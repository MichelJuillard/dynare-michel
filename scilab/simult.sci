function [y_]=simult(dr,replic,iorder,istoch,seed)
[nargout,nargin] = argn(0)
// Copyright (C) 2001 Michel Juillard
// 
// 
 
sd = 0;
if nargin>4 then
  sd = 1;
end
 
xlen = xkmax_+xkmin_+1;
klen = ykmin_+ykmax_+1;
iyv = iy_';
iyv = iyv(:);
iyr0 = find(iyv)';
it_ = ykmin_+1;
 
if exo_nbr==0 then
  exe_ = [];
end
 
if replic>1 then
  fname = fname_+'_simul';
  [fh,%v1] = mopen(fname,'wb')
  if %v1<0 then dyn_disp(%v1); end
end
 
for i = 1:replic
  if sd then
    rand('seed',seed+i-1,'n');
  end
  ex_ = rand(xkmin_+xkmax_+iter_,exo_nbr,'n')*chol(Sigma_e_);
  y_ = simult_(dr('ys'),dr,ex_,iorder,iter_);
  if replic>1 then
    mput(y_(:,ykmin_+1:$),'l',fh);
  end
end
 
if replic>1 then
  mclose(fh);
end
 
 
// 02/20/01 MJ replaced ys by dr.ys
// 02/22/01 MJ removed commented out lines
//             removed useless temps
//             stderr_ replaced by Sigma_e_
// 02/28/01 MJ changed expression for Sigma_e_
 
