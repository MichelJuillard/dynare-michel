% Copyright (C) 2001 Michel Juillard
%
function dr = olr1(ys,algo,olr_inst,bet,obj_var,W)

global M_ options_ oo_
global  it_ means_ stderrs_

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;


if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('OLR: Error in model specification: some variables don"t appear as current') ;
end

if M_.maximum_lead == 0
  error ('Backward or static model: no point in using OLR') ;
end

if xlen > 1
  error (['OLR: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end

% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = oo_.exo_steady_state';
fh = str2func([M_.fname '_static']);
if max(abs(feval(fh,ys))) > options_.dynatol
  [dr.ys, check] = dynare_solve([M_.fname '_static'],ys);
  if check
    error('OLR: convergence problem in DYNARE_SOLVE')
  end
else
  dr.ys = ys;
end
dr = olr2(dr,olr_inst,bet,obj_var,W);
oo_.exo_simul = tempex;
tempex = [];

% 04/13/03 MJ



