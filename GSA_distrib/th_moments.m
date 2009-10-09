% Copyright (C) 2001 Michel Juillard
%
function [vdec, corr, autocorr, z, zz] = th_moments(dr,var_list)
  global M_ oo_ options_
  
  nvar = size(var_list,1);
  if nvar == 0
    nvar = length(dr.order_var);
    ivar = [1:nvar]';
  else
    ivar=zeros(nvar,1);
    for i=1:nvar
      i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
      	error (['One of the variable specified does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end
  
  [gamma_y,ivar] = th_autocovariances(dr,ivar,M_, options_);
  m = dr.ys(ivar);

  
  i1 = find(abs(diag(gamma_y{1})) > 1e-12);
  s2 = diag(gamma_y{1});
  sd = sqrt(s2);

  
  z = [ m sd s2 ];
  mean = m;
  var = gamma_y{1};
  

%'THEORETICAL MOMENTS';
%'MEAN','STD. DEV.','VARIANCE');
z;

%'VARIANCE DECOMPOSITION (in percent)';
if M_.exo_nbr>1,
vdec = 100*gamma_y{options_.ar+2}(i1,:);
else
vdec = 100*ones(size(gamma_y{1}(i1,1)));
end  
%'MATRIX OF CORRELATIONS';
    corr = gamma_y{1}(i1,i1)./(sd(i1)*sd(i1)');
    corr = corr-diag(diag(corr))+diag(diag(gamma_y{1}(i1,i1)));
  
  if options_.ar > 0
%'COEFFICIENTS OF AUTOCORRELATION';
    for i=1:options_.ar
      autocorr{i} = gamma_y{i+1};
      zz(:,i) = diag(gamma_y{i+1}(i1,i1));
    end
  end
  

  