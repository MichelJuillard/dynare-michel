% Copyright (C) 2001 Michel Juillard
%
function steady_()

  global M_ oo_ it_ options_

  if exist([M_.fname '_steadystate'])
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],...
				     oo_.steady_state,...
				     [oo_.exo_steady_state; ...
		                      oo_.exo_det_steady_state]);
  else
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
				     oo_.steady_state,...
				     options_.jacobian_flag, ...	    
			             [oo_.exo_steady_state; ...
		                      oo_.exo_det_steady_state]);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end


