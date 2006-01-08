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

  if M_.exo_det_nbr > 0
    tempdet = oo_.exo_det_simul ;
    oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+1,1) ;
  end

  if exist([M_.fname '_steadystate'])
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],x);
  else
    [oo_.staedy_state,check] = dynare_solve([M_.fname '_static'],x);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end

  if M_.exo_det_nbr > 0
    oo_.exo_det_simul = tempdet;
  end
  oo_.exo_simul = temp ;
