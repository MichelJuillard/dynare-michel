% Copyright (C) 2001 Michel Juillard
%
function resid(options_.periods)
  global M_ options_ oo_ it_ endval_ z
  
  oo_.exo_simul = ones(M_.maximum_lag+M_.maximum_lead+options_.periods,1)*oo_.exo_steady_state';
  n = size(M_.lead_lag_incidence,2);
%  if ~ options_.initval_file | size(oo_.y_simul,2) ~= options_.periods+M_.maximum_lag+M_.maximum_lead
  if ~ options_.initval_file 
    if size(oo_.steady_state,1) == 1 & oo_.steady_state == 0
      oo_.steady_state = zeros(size(oo_.steady_state,1),1) ;
    end
    oo_.y_simul = oo_.steady_state*ones(1,options_.periods+M_.maximum_lag+M_.maximum_lead) ;
    if endval_ == 1
      oo_.y_simul(:,1:M_.maximum_lag) = ys0_*ones(1,M_.maximum_lag) ;
    end
  end

  i = M_.lead_lag_incidence';
  iyr0 = find(i(:));

  y =oo_.y_simul(:);
  z = zeros(n,options_.periods);
  fh = str2func([M_.fname '_static']);
  for it_=M_.maximum_lag+1:options_.periods+M_.maximum_lag
    z(:,it_-M_.maximum_lag) = feval(fh,y(iyr0));
    iyr0 = iyr0 + n;
  end

  disp([[1:options_.periods]' z']); 






