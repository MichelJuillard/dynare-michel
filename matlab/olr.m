% Copyright (C) 2001 Michel Juillard
%
function olr(var_list,olr_inst,obj_var,W)
  global M_ options_ oo_ 
    
  options_.order = 1;
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
  options_ = set_default_option(options_,'olr_beta',1);
  options_ = set_default_option(options_,'periods',1);
  if options_.simul & ~isempty(options_.periods) & options_.periods == 0
    options_.periods = options_.periods;
  end

  options_.periods = max(options_.periods,1);
  
  make_ex_;
  
  oo_.dr=olr1(oo_.steady_state,options_.dr_algo,olr_inst,options_.olr_beta,obj_var,W);

  disp(' ')
  disp('OPTIMAL LINEAR REGULATOR')
  disp(' ')
  disp('MODEL SUMMARY')
  disp(' ')
  disp(['  Number of variables:         ' int2str(M_.endo_nbr-size(olr_inst,1))])
  disp(['  Number of stochastic shocks: ' int2str(M_.exo_nbr)])
  disp(['  Number of instruments        ' int2str(size(olr_inst,1))])
%  disp(['  Number of state variables:   ' ...
%	int2str(length(find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1)))])
%  disp(['  Number of jumpers:           ' ...
%	int2str(length(find(oo_.dr.kstate(:,2) == M_.maximum_lag+2)))])
  disp(['  Number of static variables:  ' int2str(oo_.dr.nstatic)])
  my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
  labels = deblank(M_.exo_name);
  headers = strvcat('Variables',labels);
  lh = size(labels,2)+2;
  table(my_title,headers,labels,M_.Sigma_e,lh,10,6);
  disp(' ')
  disp_dr(oo_.dr,options_.order,var_list);
  if options_.order == 1 & options_.simul == 0 & options_.nomoments == 0
    disp_th_moments(oo_.dr,var_list);
  elseif options_.simul == 1 | options_.nomoments == 0
    if options_.periods == 0
      error('OLR error: number of periods for the simulation isn''t specified')
    end
    if options_.periods < options_.drop
      disp(['OLR error: The horizon of simulation is shorter' ...
	    ' than the number of observations to be DROPed'])
      return
    end
    
    oo_.y_simul = simult(repmat(oo_.dr.ys,1,M_.maximum_lag),oo_.dr);
    dyn2vec;
    if options_.nomoments == 0
      disp_moments(oo_.y_simul,var_list);
    end
  end
  
  n = size(var_list,1);
  if n == 0
    n = length(oo_.dr.order_var);
    ivar = [1:n]';
    var_list = M_.endo_names;
  else
    ivar=zeros(n,1);
    for i=1:n
      i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
      	error (['One of the specified variables does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end

  if n < 13 & options_.irf > 0
    if n == 1
      nr = 1;
      nc = 1;
    elseif n == 2
      nr = 1;
      nc = 2;
    elseif n <= 4
      nr = 2;
      nc = 2;
    elseif n <= 6
      nr = 2;
      nc = 3;
    elseif n <= 9
      nr = 3;
      nc = 3;
    elseif n <= 12
      nr = 3;
      nc = 4;
    end
    olditer = options_.periods;
    if options_.order == 1
      options_.replic = 1;
    else
      if options_.replic == 0
	options_.replic = 50;
      end
    end
    for i = 1:M_.exo_nbr
      figure('Name',['Shock to ' M_.exo_name(i,:)]);
      y=irf(oo_.dr,M_.exo_name(i,:),sqrt(M_.Sigma_e(i,i)), options_.irf, options_.drop, options_.replic, options_.order);
      m = 1;
      for j = 1:n
	if max(y(ivar(j),:))-min(y(ivar(j),:)) > 1e-10
	  subplot(nr,nc,m);
	  plot([y(ivar(j),:)']);
	  title(var_list(j,:),'Interpreter','none');
	  assignin('base',[deblank(var_list(j,:)) '_' deblank(M_.exo_name(i,:))],y(ivar(j),:)');
	  m = m + 1;
	end
      end
    
    end
    options_.periods = olditer;
  end
      
% 05/21/03 MJ