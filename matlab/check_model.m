function check_model()
  global M_
  
  xlen = M_.maximum_exo_lag+M_.maximum_exo_lead + 1;
  if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('RESOL: Error in model specification: some variables don"t appear as current') ;
end

if xlen > 1
  error (['RESOL: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end

if (M_.exo_det_nbr > 0) & (M_.maximum_lag > 1 | M_.maximum_lead > 1)
  error(['Exogenous deterministic variables are currently only allowed in' ...
	 ' models with leads and lags on only one period'])
end

