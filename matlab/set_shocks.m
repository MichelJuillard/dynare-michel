% Copyright (C) 2003 Michel Juillard
%
% set
function set_shocks(flag,k,ivar,values)
  global oo_
  
  n = size(oo_.exo_simul,1);
  if k(end) > n
    oo_.exo_simul = [oo_.exo_simul; ones(k(end)-n,1)*oo_.exo_steady_state'];
  end
  
  if flag == 0
    oo_.exo_simul(k,ivar) = ones(length(k),1).*values;
  else
    oo_.exo_simul(k,ivar) = oo_.exo_simul(k,ivar).*values;
  end

  % 05/29/03 MJ