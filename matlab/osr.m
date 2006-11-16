% Copyright (C) 2001 Michel Juillard
%
function osr(var_list,params,i_var,W)
  global M_ options_ oo_  

  options_.order = 1;
  options_ = set_default_option(options_,'linear',0);
  options_ = set_default_option(options_,'ar',5);
  options_ = set_default_option(options_,'irf',40);
  options_ = set_default_option(options_,'dr_algo',0);
  options_ = set_default_option(options_,'simul_algo',0);
  options_ = set_default_option(options_,'drop',100);
  options_ = set_default_option(options_,'replic',1);
  options_ = set_default_option(options_,'nomoments',0);
  options_ = set_default_option(options_,'nocorr',0);
  options_ = set_default_option(options_,'simul_seed',[]);
  options_ = set_default_option(options_,'hp_filter',0);
  options_ = set_default_option(options_,'hp_ngrid',512);
  options_ = set_default_option(options_,'simul',0);
  options_ = set_default_option(options_,'periods',1);
  if options_.simul & ~isempty(options_.periods) & options_.periods == 0
    options_.periods = options_.periods;
  end

  options_.periods = max(options_.periods,1);
  options_.periods = options_.periods;
  
  make_ex_;

  np = size(params,1);
  i_params = zeros(np,1);
  for i=1:np
    i_params(i) = strmatch(deblank(params(i,:)),M_.param_names,'exact');
  end
    
  disp(' ')
  disp('OPTIMAL SIMPLE RULE')
  disp(' ')
  osr1(i_params,i_var,W);

  stoch_simul(var_list);