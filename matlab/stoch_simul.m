% Copyright (C) 2001 Michel Juillard
%
function stoch_simul(var_list)
  global M_ options_ oo_   
    
  options_ = set_default_option(options_,'order',2);
  options_ = set_default_option(options_,'linear',0);
  if options_.linear
    options_.order = 1;
  end
  options_ = set_default_option(options_,'ar',5);
  options_ = set_default_option(options_,'irf',40);
  options_ = set_default_option(options_,'dr_algo',0);
  options_ = set_default_option(options_,'simul_algo',0);
  options_ = set_default_option(options_,'drop',100);
  if options_.order == 1
    options_.replic = 1;
  else
    options_ = set_default_option(options_,'replic',50);
  end
  options_ = set_default_option(options_,'nomoments',0);
  options_ = set_default_option(options_,'nocorr',0);
  options_ = set_default_option(options_,'simul_seed',[]);
  options_ = set_default_option(options_,'hp_filter',0);
  options_ = set_default_option(options_,'hp_ngrid',512);
  options_ = set_default_option(options_,'simul',0);
  options_ = set_default_option(options_,'periods',0);
  if options_.simul & ~isempty(options_.periods) & options_.periods == 0
    options_.periods = options_.periods;
  end
  options_.periods = max(options_.periods,1);

  if M_.exo_nbr > 0
    oo_.exo_simul= repmat(oo_.exo_steady_state',options_.periods + M_.maximum_lag + M_.maximum_lead,1) ;
  end

  oo_.dr=resol(oo_.steady_state,options_.dr_algo,options_.linear,options_.order);

  disp(' ')
  disp('MODEL SUMMARY')
  disp(' ')
  disp(['  Number of variables:         ' int2str(M_.endo_nbr)])
  disp(['  Number of stochastic shocks: ' int2str(M_.exo_nbr)])
  disp(['  Number of state variables:   ' ...
	int2str(length(find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1)))])
  disp(['  Number of jumpers:           ' ...
	int2str(length(find(oo_.dr.kstate(:,2) == M_.maximum_lag+2)))])
  disp(['  Number of static variables:  ' int2str(oo_.dr.nstatic)])
  my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
  labels = deblank(M_.exo_names);
  headers = strvcat('Variables',labels);
  lh = size(labels,2)+2;
  table(my_title,headers,labels,M_.Sigma_e,lh,10,6);
  disp(' ')
  disp_dr(oo_.dr,options_.order,var_list);
  
  if options_.simul == 0 & options_.nomoments == 0
    disp_th_moments(oo_.dr,var_list);
  elseif options_.simul == 1
    if options_.periods == 0
      error('STOCH_SIMUL error: number of periods for the simulation isn''t specified')
    end
    if options_.periods < options_.drop
      disp(['STOCH_SIMUL error: The horizon of simulation is shorter' ...
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
    n = M_.endo_nbr;
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
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    cs = chol(SS)';
    tit(M_.exo_names_orig_ord,:) = M_.exo_names;
    for i = 1:M_.exo_nbr
      if SS(i,i) > 1e-13
	figure('Name',['Orthogonalized shock to ' tit(i,:)]);
	y=irf(oo_.dr,cs(M_.exo_names_orig_ord,i), options_.irf, options_.drop, ...
	      options_.replic, options_.order);
	m = 1;
	for j = 1:n
	  if max(y(ivar(j),:)) - min(y(ivar(j),:)) > 1e-10 
	    subplot(nr,nc,m);
	    plot([y(ivar(j),:)']);
	    title(var_list(j,:),'Interpreter','none');
	    assignin('base',[deblank(var_list(j,:)) '_' deblank(tit(i,:))],y(ivar(j),:)');
	    m = m+1;
	  end
	end
      end
    end
    options_.periods = olditer;
  end
      
% 01/10/01 FC oo_.dr and oo_.y_simul made global
% 02/20/01 MJ oo_.steady_state removed from calling sequence for simult (all in oo_.dr)
% 02/23/01 MJ added dyn2vec()
% 06/24/01 MJ steady -> steady_
% 09/24/01 MJ oo_.dr made global
% 08/28/02 MJ added var_list
% 10/09/02 MJ no simulation and theoretical moments for order 1 
% 10/14/02 MJ added plot of IRFs
% 10/30/02 MJ options_ are now a structure
% 01/01/03 MJ added dr_algo
% 01/09/03 MJ set default values for options_ (correct absence of autocorr
%             when order == 1)
% 01/12/03 MJ removed call to steady_ as already checked in resol()
% 01/31/03 MJ make IRF global with varname_shockname
% 02/09/03 MJ oo_.steady_state reset with value declared in initval after computations
% 02/18/03 MJ removed above change. oo_.steady_state shouldn't be affected by
%             computations in this function
%             new option SIMUL computes a stochastic simulation and save
%             results in oo_.y_simul and via dyn2vec
% 04/03/03 MJ corrected bug for simulation with M_.maximum_lag > 1
% 05/20/03 MJ eliminates exogenous shocks with 0 variance
% 05/20/03 MJ don't plot IRF if variation < 1e-10
% 11/14/03 MJ corrected bug on number of replications for IRF when
%             order=2
% 11/22/03 MJ replaced IRFs by orthogonalized IRFs
