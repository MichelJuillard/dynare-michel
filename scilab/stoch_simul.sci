function []=stoch_simul(options,var_list)
// Copyright (C) 2001 Michel Juillard
// 
global y_ dr_ ys_
 
if options.linear
  options.order = 1;
end
if isempty(options.ar)
  options.ar = 5;
end
if isempty(options.irf)
  options.irf = 40;
end
dr_algo = options.dr_algo;
simul_algo = options.simul_algo;
drop = options.drop;
%linear = options.linear;
replic = options.replic;
irf_length = options.irf;
order = options.order;
if %linear == 1 then
  order = 1;
else
  order = options.order;
end
if iter_ < drop then
  dyn_disp('STOCH_SIMUL error: The horizon of simulation is shorter than the number'+' of observations to be DROPed');
  return
end

temps = ys_;
dr_ = resol(ys_,dr_algo,%linear,order);
disp_dr(dr_,order,var_list);
if options.nomoments == 0
  if order == 1 then
    disp_th_moments(dr_,var_list);
  else
    y_ = simult(dr_, 1, order, 1);
    ys_ = temps;
    dyn2vec();
    disp_moments(y_,drop,var_list);
  end
end

n = size(var_list,1);
if n == 0 then
n = endo_nbr;
ivar = [1:n]';
var_list = lgy_;
else
  ivar=zeros(n,1);
  for i=1:n
    i_tmp = grep_exact(lgy_,var_list(i,:));
    if isempty(i_tmp) then
      error (['One of the specified variables does not exist']) ;
    else
      ivar(i) = i_tmp;
    end
  end
end

if n < 16 & options.irf > 0 then
  if n == 1 then
    nr = 1;
    nc = 1;
  elseif n == 2 then
    nr = 1;
    nc = 2;
  elseif n <= 4 then
    nr = 2;
    nc = 2;
  elseif n <= 6 then
    nr = 2;
    nc = 3;
  elseif n <= 9 then
    nr = 3;
    nc = 3;
  elseif n <= 12 then
    nr = 3;
    nc = 4;
  elseif n <= 16 then
    nr = 4;
    nc = 4;
  end
  olditer = iter_;
  for i = 1:exo_nbr
    xset('window',i-1);
    xname('Shock to '+lgx_(i));
    if order == 1
      replic = 1;
    elseif replic == 0; 
      replic = 100;
    end
    y=irf(dr_,lgx_(i,:),sqrt(Sigma_e_(i,i)), irf_length, drop, replic, order);
    for j = 1:n
      subplot(nr,nc,j);
      plot2d([y(ivar(j),:)']);
      xtitle(var_list(j,:));
    end
  end
  iter_ = olditer;
end
ys_ = temps;
 
// 01/10/01 FC dr_ and y_ made global
// 02/20/01 MJ ys_ removed from calling sequence for simult (all in dr_)
// 02/23/01 MJ added dyn2vec()
// 06/24/01 MJ steady -> steady_
// 09/24/01 MJ dr_ made global
// 01/14/03 MJ changed options passing and defautls.
//             replaced algo by dr_algo and simul_algo
 
 
