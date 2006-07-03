function []=disp_moments(y,drop,var_list)
// Copyright (C) 2001 Michel Juillard
// 

nvar = size(var_list,1);
if nvar == 0
  nvar = endo_nbr;
  ivar = [1:nvar]';
else
  ivar=zeros(nvar,1);
  for i=1:nvar
    i_tmp = grep_exact(lgy_,var_list(i,:));
    if isempty(i_tmp)
      error (['One of the variable specified does not exist']) ;
    else
      ivar(i) = i_tmp;
    end
  end
end
y = y(ivar,drop+1:$)';

m = mean(y,1);
y = y-ones(size(y,1),1).*.m;
s2 = mean(y .* y,1);
z = [m',s2',(mean(y.^3,1) ./ (s2.^1.5))',(mean(y.^4,1) ./ (s2 .* s2)-3)'];
 
dyn_disp('MOMENTS OF SIMULATED VARIABLES');
dyn_disp('VARIABLE     MEAN       VARIANCE   SKEWNESS  KURTOSIS');
for i = 1:nvar
  mprintf('%10s %10.6f %10.6f %10.6f %10.6f\n',lgy_(i,:),z(i,:))
  mfprintf(fh_log,'%10s %10.6f %10.6f %10.6f %10.6f\n',lgy_(i,:),z(i,:))
end
 
