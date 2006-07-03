function []=disp_th_moments(dr,var_list)
// Copyright (C) 2001 Michel Juillard
// 
global('endo_nbr','lgy_');
 
nvar = size(var_list,1);
order = dr('order_var');
if nvar==0 then
  nvar = endo_nbr;
  ivar = (1:nvar)';
else
  ivar = zeros(nvar,1);
  for i = 1:nvar
     
    //!! Unknown function strmatch ,the original calling sequence is used
    i_tmp = grep_exact(lgy_,var_list(i,:));
    if i_tmp==[] then
      error('One of the variable specified does not exist');
    else
      ivar(i,1) = i_tmp(:);
    end
  end
end
 
m = dr('ys')(ivar);
 
 
//!! Unknown function th_autocovariances ,the original calling sequence is used
Gamma_y = th_autocovariances(dr,ivar);
 
s2 = diag(Gamma_y);
z = [m,sqrt(s2),s2];
 
dyn_disp('THEORETICAL MOMENTS');
dyn_disp('VARIABLE           MEAN             STD. DEV.           VARIANCE');
for i = 1:size(ivar,1)
  dyn_disp(sprintf('%16s %16.6f %16.6f %16.6f\n',lgy_(ivar(i),:),z(i,:)));
end
 
// 10/09/02 MJ 
// 10/18/02 MJ added th_autocovariances() and provided for lags on several
// periods
