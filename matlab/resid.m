% Copyright (C) 2001 Michel Juillard
%
function resid(period)
  global M_ options_ oo_ it_ endval_ z
  
  oo_.exo_simul = ones(M_.maximum_lag+M_.maximum_lead+period,1)*oo_.exo_steady_state';
  n = size(M_.lead_lag_incidence,2);
%  if ~ options_.initval_file | size(oo_.endo_simul,2) ~= period+M_.maximum_lag+M_.maximum_lead
  if ~ options_.initval_file 
    if size(oo_.steady_state,1) == 1 & oo_.steady_state == 0
      oo_.steady_state = zeros(size(oo_.steady_state,1),1) ;
    end
    oo_.endo_simul = oo_.steady_state*ones(1,period+M_.maximum_lag+M_.maximum_lead) ;
    if endval_ == 1
      oo_.endo_simul(:,1:M_.maximum_lag) = ys0_*ones(1,M_.maximum_lag) ;
    end
  end

  i = M_.lead_lag_incidence';
  iyr0 = find(i(:));

  y =oo_.endo_simul(:);
  z = zeros(n,period);
  fh = str2func([M_.fname '_dynamic']);
  for it_=M_.maximum_lag+1:period+M_.maximum_lag
    z(:,it_-M_.maximum_lag) = feval(fh,y(iyr0),oo_.exo_simul);
    iyr0 = iyr0 + n;
  end

  % disp([[1:period]' z']); 

  for i = 1:4
    disp(' ')
  end
  
  for i=1:length(z)
    if abs(z(i)) < 10^(-12)
      z(i) = 0;
    end
    disp(['Residual for equation number ' int2str(i) ' is equal to ' num2str(z(i))])
  end  
  
  for i = 1:2
    disp(' ')
  end




