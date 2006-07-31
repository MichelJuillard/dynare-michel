function osr1(i_params,weights)
  global M_ oo_ options_ it_

  klen = M_.maximum_lag + M_.maximum_lead + 1;
  iyv = M_.lead_lag_incidence';
  iyv = iyv(:);
  iyr0 = find(iyv) ;
  it_ = M_.maximum_lag + 1 ;


  if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
  end

  if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
    error ('OSR: Error in model specification: some variables don"t appear as current') ;
  end

  if M_.maximum_lead == 0
    error ('Backward or static model: no point in using OLR') ;
  end

  exe =zeros(M_.exo_nbr,1);
  % check if ys is steady state
  fh = str2func([M_.fname '_static']);
  if max(abs(feval(fh,oo_.steady_state,exe))) > options_.dynatol
    [oo_.dr.ys, check] = dynare_solve([M_.fname '_static'],oo_.steady_state,exe);
    if check
      error('OLR: convergence problem in DYNARE_SOLVE')
    end
  else
    dr.ys = oo_.steady_state;
    oo_.dr = dr;
  end

  
  np = size(i_params,1);
  t0 = M_.params(i_params);
  
  options = optimset('fminunc');
  options = optimset('display','iter');
  [p,f]=fminunc(@osr_obj,t0,options,i_params,weights);


  
  disp('')
  disp('OPTIMAL VALUE OF THE PARAMETERS:')
  disp('')
  for i=1:np
    disp(sprintf('%16s %16.6g\n',M_.param_names(i_params(i),:),p(i)))
  end
  disp(sprintf('Objective function : %16.6g\n',f));
  disp(' ')
  oo_.dr=resol(oo_.steady_state,0);

  % 05/10/03 MJ modified to work with osr.m and give full report