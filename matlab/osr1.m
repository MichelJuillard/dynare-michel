function oo_.dr=osr1(params,weights)
  global M_ oo_ options_ it_

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
  fh = str2func([M_.fname '_static']);
  if max(abs(feval(fh,oo_.steady_state))) > options_.dynatol
    [oo_.dr.ys, check] = dynare_solve([M_.fname '_static'],oo_.steady_state);
    if check
      error('OLR: convergence problem in DYNARE_SOLVE')
    end
  else
    oo_.dr.ys = oo_.steady_state;
  end

  
  np = size(params,1);
  t0 = zeros(np,1);
  for i=1:np
    t0(i)=evalin('base',[params(i,:) ';']);
  end
  
  [p,f]=fminsearch(@osr_obj,t0,[],params,weights);

  disp('')
  disp('OPTIMAL VALUE OF THE PARAMETERS:')
  disp('')
  for i=1:np
    disp(sprintf('%16s %16.6g\n',params(i,:),p(i)))
  end
  disp(sprintf('Objective function : %16.6g\n',f));
  disp(' ')
  oo_.dr=resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);

  % 05/10/03 MJ modified to work with osr.m and give full report