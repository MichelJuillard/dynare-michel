% Copyright (C) 2001 Michel Juillard
%
function steady_()

  global M_ oo_ it_


  x = oo_.steady_state ;
  xlen = M_.maximum_lag + M_.maximum_lead + 1 ;
  nn = size(M_.lead_lag_incidence,2) ;
  it_ = M_.maximum_lag+1 ;
  temp = oo_.exo_simul ;
  oo_.exo_simul = ones(xlen,1)*transpose(oo_.exo_steady_state);

  [oo_.steady_state,cheik] = dynare_solve([M_.fname '_static'],x,oo_.exo_simul);

  if cheik ~= 0
    error('STEADY: convergence problems')
  end

  oo_.exo_simul = temp ;

% 06/24/01 MJ: steady_ no results printer; steady with printed results
% 07/31/03 MJ: in case of convergence problem steady stops with an error
% 01/22/05 SA: check --> cheik (to avoid confusion with function check.m)
