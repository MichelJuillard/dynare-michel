function set_shocks(flag,k,ivar,values)

% function set_shocks(flag,k,ivar,values)
% writes a deterministic shock into exo_simul or exo_det_simul
%
% INPUTS
%    flag=0:    replaces exo_simul  
%    flag=1:    multiplicative shock into exo_simul
%    flag=2:    replaces exo_det_simul
%    flag=3:    multipliczative shock into exo_det_simul
%    k:         period
%    ivar:      indice of exogenous variables
%    values:    shock values
%
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  global oo_ M_
  
  k = k + M_.maximum_lag;
  n1 = size(oo_.exo_simul,1);
  n2 = size(oo_.exo_det_simul,1);
  if k(end) > n1 & flag <= 1
    oo_.exo_simul = [oo_.exo_simul; repmat(oo_.exo_steady_state',k(end)-n1,1)];
  elseif k(end) > n2 & flag > 1
    oo_.exo_det_simul = [oo_.exo_det_simul; repmat(oo_.exo_det_steady_state',k(end)-n2,1)];
  end
  
  switch flag
   case 0
    oo_.exo_simul(k,ivar) = repmat(values,length(k),1);
   case 1
    oo_.exo_simul(k,ivar) = oo_.exo_simul(k,ivar).*values;
   case 2
    oo_.exo_det_simul(k,ivar) = repmat(values,length(k),1);
   case 3
    oo_.exo_det_simul(k,ivar) = oo_.exo_det_simul(k,ivar).*values;
  end

