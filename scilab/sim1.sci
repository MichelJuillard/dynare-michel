function []=sim1()
// Copyright (C) 2001 Michel Juillard
// 
 
global it_ jacobia_ d1_ y_ c

ny = size(y_,1);
nyp = nnz(iy_(1,:));
nyf = nnz(iy_(3,:));
nrs = ny+nyp+nyf+1;
nrc = nyf+1;
iyf = find(iy_(3,:)>0);
iyp = find(iy_(1,:)>0);
isp = 1:nyp;
is = nyp+1:ny+nyp;
isf = iyf+nyp;
isf1 = nyp+ny+1:nyf+nyp+ny+1;
stop = 0;
 
dyn_disp('-----------------------------------------------------');
dyn_disp('MODEL SIMULATION :');
dyn_disp('');

if start_simul==[] then
  it_init = 2;
else
  it_init = start_simul;
end
 
h1 = 0;
for iter = 1:maxit_
  timer();
   
  if ct_==0 then
    c = zeros(ny*iter_,nrc);
  else
    c = zeros(ny*(iter_+1),nrc);
  end
   
  it_ = it_init;
  z = [y_(iyp,it_-1);y_(:,it_);y_(iyf,it_+1)];
  jacob(fname_+'_ff',z);
  jacobia_ = [jacobia_,-d1_];
  ic = 1:ny;
  icp = iyp;
  c(ic,:) = jacobia_(:,is)\jacobia_(:,isf1);
  for it_ = it_init+1:iter_+1
    z = [y_(iyp,it_-1);y_(:,it_);y_(iyf,it_+1)];
    jacob(fname_+'_ff',z);
    jacobia_ = [jacobia_,-d1_];
    jacobia_(:,[isf,nrs]) = jacobia_(:,[isf,nrs])-jacobia_(:,isp)*c(icp,:);
    ic = ic+ny;
    icp = icp+ny;
    c(ic,:) = jacobia_(:,is)\jacobia_(:,isf1);
  end
   
  if ct_==1 then
    s = eye(ny,ny);
    s(:,isf) = s(:,isf)+c(ic,1:nyf);
    ic = ic+ny;
    c(ic,nrc) = s\c(:,nrc);
    c = bksup1(ny,nrc);
    c = matrix(c,ny,iter_+1);
    y_(:,it_init:iter_+2) = y_(:,it_init:iter_+2)+slowc_*c;
  else
    c = bksup1(ny,nrc);
    c = matrix(c,ny,iter_);
    y_(:,it_init:iter_+1) = y_(:,it_init:iter_+1)+slowc_*c;
  end
   
  err = max(max(abs(c ./ scalv_'),'r'));
  dyn_disp(string(iter)+' - err = '+string(err));
  
  h2 = timer();
  dyn_disp(' Time of iteration  :'+string(h2));
  h1 = h1 + h2;
   
  if err < dynatol_ then
    stop = 1;
    
    dyn_disp(''); 
    dyn_disp(' Total time of simulation  :'+string(h1));
    dyn_disp('');
    dyn_disp(' Convergency obtained.');
    dyn_disp('');
    break
     ;
     
  end
end
 
if ~stop then
  dyn_disp('');
  dyn_disp(' Total time of simulation  :'+string(h1));
  dyn_disp('');
  dyn_disp('WARNING : maximum number of iterations is reached (modify maxit_).');
  dyn_disp('');
end
dyn_disp('-----------------------------------------------------');
return
 
 
// 08/24/01 MJ added start_simul
// 09/26/01 MJ translated to Scilab
