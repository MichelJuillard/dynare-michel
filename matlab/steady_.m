% Copyright (C) 2001 Michel Juillard
%
function steady_()

  global M_ oo_ it_


  x = oo_.steady_state ;
  xlen = M_.maximum_lag + M_.maximum_lead + 1 ;
  nn = size(M_.lead_lag_incidence,2) ;
  it_ = M_.maximum_lag+1 ;
  x = repmat(oo_.exo_steady_state',xlen,1);

  if M_.exo_det_nbr > 0
    x = [x, repmat(oo_.exo_det_steady_state',M_.maximum_lag+1,1)] ;
  end

  if exist([M_.fname '_steadystate'])
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],oo_.steady_state,x);
  else
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],oo_.steady_state,1,x);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end


