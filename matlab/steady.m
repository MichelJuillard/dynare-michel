% Copyright (C) 2001 Michel Juillard
%
function steady(linear)

  global M_ oo_ options_ ys0_ 

  options_ = set_default_option(options_,'jacobian_flag',1);
  options_ = set_default_option(options_,'steadystate_flag',0);
  if exist([M_.fname '_steadystate.m'])
    options_.steadystate_flag = 1;
  end 
    
  steady_;
  
  disp(' ')
  disp('STEADY-STATE RESULTS:')
  disp(' ')
  endo_names = M_.endo_names;
  steady_state = oo_.steady_state;
  for i=1:size(oo_.steady_state,1)
    disp(sprintf('%s \t\t %g',endo_names(i,:),steady_state(i)));
  end
  
  if isempty(ys0_)
    oo_.y_simul(:,1:M_.maximum_lag) = oo_.steady_state * ones(1,M_.maximum_lag);
  else
    options_ =set_default_option(options_,'periods',1);
    oo_.y_simul(:,M_.maximum_lag+1:M_.maximum_lag+options_.periods+M_.maximum_lead) = oo_.steady_state * ones(1,options_.periods+M_.maximum_lead);
  end
  
% 06/24/01 MJ steady print results; steady_ doesn't
% 09/22/01 FC corrected lgy(i,:)
% 05/29/03 MJ sets initial values of oo_.y_simul
