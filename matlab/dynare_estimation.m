function dynare_estimation(var_list_)

global M_ options_ oo_ estim_params_ 
global bayestopt_ dsge_prior_weight

% temporary fix until M_.H is initialized by the parser
M_.H = [];

options_.varlist = var_list_;
options_.lgyidx2varobs = zeros(size(M_.endo_names,1),1);
for i = 1:size(M_.endo_names,1)
  tmp = strmatch(deblank(M_.endo_names(i,:)),options_.varobs,'exact');
  if ~isempty(tmp)
    options_.lgyidx2varobs(i,1) = tmp;
  end
end  

options_ = set_default_option(options_,'first_obs',1);
options_ = set_default_option(options_,'prefilter',0);
options_ = set_default_option(options_,'presample',0);
options_ = set_default_option(options_,'lik_algo',1);
options_ = set_default_option(options_,'lik_init',1);
options_ = set_default_option(options_,'nograph',0);
options_ = set_default_option(options_,'mh_conf_sig',0.90);
options_ = set_default_option(options_,'mh_replic',20000);
options_ = set_default_option(options_,'mh_drop',0.5);
options_ = set_default_option(options_,'mh_jscale',0.2);
options_ = set_default_option(options_,'mh_init_scale',2*options_.mh_jscale);
options_ = set_default_option(options_,'mode_file','');
options_ = set_default_option(options_,'mode_compute',4);
options_ = set_default_option(options_,'mode_check',0);
options_ = set_default_option(options_,'prior_trunc',1e-10);
options_ = set_default_option(options_,'mh_mode',1); 	
options_ = set_default_option(options_,'mh_nblck',2);	
options_ = set_default_option(options_,'load_mh_file',0);
options_ = set_default_option(options_,'nodiagnostic',0);
options_ = set_default_option(options_,'loglinear',0);
options_ = set_default_option(options_,'unit_root_vars',[]);
options_ = set_default_option(options_,'XTick',[]);
options_ = set_default_option(options_,'XTickLabel',[]);
options_ = set_default_option(options_,'bayesian_irf',0);
options_ = set_default_option(options_,'bayesian_th_moments',0);
options_ = set_default_option(options_,'TeX',0);
options_ = set_default_option(options_,'irf',40);
options_ = set_default_option(options_,'relative_irf',0);
options_ = set_default_option(options_,'order',1);
options_ = set_default_option(options_,'ar',5);
options_ = set_default_option(options_,'dr_algo',0);
options_ = set_default_option(options_,'linear',0);
options_ = set_default_option(options_,'drop',0);
options_ = set_default_option(options_,'replic',1);
options_ = set_default_option(options_,'hp_filter',0);
options_ = set_default_option(options_,'smoother',0);
options_ = set_default_option(options_,'moments_varendo',0);
options_ = set_default_option(options_,'filtered_vars',0);
options_ = set_default_option(options_,'kalman_algo',1);
options_ = set_default_option(options_,'kalman_tol',10^(-12));
options_ = set_default_option(options_,'posterior_mode_estimation',1);
options_ = set_default_option(options_,'MaxNumberOfBytes',1e6);
options_ = set_default_option(options_,'xls_sheet','');
options_ = set_default_option(options_,'xls_range','');
options_ = set_default_option(options_,'filter_step_ahead',0);
options_ = set_default_option(options_,'diffuse_d',[]);
options_ = set_default_option(options_,'Opt6Iter',3);
options_ = set_default_option(options_,'Opt6Numb',100000);
options_ = set_default_option(options_,'steadystate_flag',0);
options_ = set_default_option(options_,'logdata',0);
options_ = set_default_option(options_,'use_mh_covariance_matrix',0);
options_ = set_default_option(options_,'noconstant',0);
options_ = set_default_option(options_,'steadystate_partial',[]);

if options_.order > 1
  options_.order = 1;
end

if options_.prefilter == 1
  options_.noconstant = 1;
end

if options_.filtered_vars ~= 0 & options_.filter_step_ahead == 0
  options_.filter_step_ahead = 1;
end
if options_.filter_step_ahead ~= 0
  options_.nk = max(options_.filter_step_ahead);
else
  options_.nk = 0;
end

%% Add something to the parser ++>
M_.dname = M_.fname; % The user should be able to choose another name
                     % for the directory...


pnames 		= ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
n_varobs 	= size(options_.varobs,1);

[xparam1,estim_params_,bayestopt_,lb,ub] = set_prior(estim_params_);

if any(bayestopt_.pshape > 0)
  if options_.mode_compute
    plot_priors
  end
else
  options_.mh_replic = 0;
end

bounds = prior_bounds(bayestopt_);
bounds(:,1)=max(bounds(:,1),lb);
bounds(:,2)=min(bounds(:,2),ub);

if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
  find(xparam1 < bounds(:,1))
  find(xparam1 > bounds(:,2))
  error('Initial parameter values are outside parameter bounds')
end
lb = bounds(:,1);
ub = bounds(:,2);
bayestopt_.lb = lb;
bayestopt_.ub = ub;

if ~isfield(options_,'trend_coeffs')
  bayestopt_.with_trend = 0;
else
  bayestopt_.with_trend = 1;
  bayestopt_.trend_coeff = {};
  trend_coeffs = options_.trend_coeffs;
  nt = length(trend_coeffs);
  for i=1:n_varobs
    if i > length(trend_coeffs)
      bayestopt_.trend_coeff{i} = '0';
    else
      bayestopt_.trend_coeff{i} = trend_coeffs{i};
    end
  end
end

bayestopt_.penalty = 1e8;	% penalty 

nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np ;
nx  = nvx+nvn+ncx+ncn+np;

dr = set_state_space([]);
nstatic = dr.nstatic;
npred = dr.npred;
nspred = dr.nspred;

if isempty(options_.varobs)
  error('ESTIMATION: VAROBS is missing')
end

%% Setting resticted state space (observed + predetermined variables)

k = [];
k1 = [];
for i=1:n_varobs
  k = [k strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:),'exact')];
  k1 = [k1 strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
end
% union of observed and state variables
k2 = union(k',[dr.nstatic+1:dr.nstatic+dr.npred]');

% set restrict_state to postion of observed + state variables
% in expanded state vector
bayestopt_.restrict_var_list = k2;
% set mf1 to positions of observed variables in restricted state vector
% for likelihood computation
[junk,bayestopt_.mf1] = ismember(k,k2); 
% set mf2 to positions of observed variables in expanded state vector
% for filtering and smoothing
bayestopt_.mf2 	= k;
bayestopt_.mfys = k1;

[junk,ic] = intersect(k2,nstatic+(1:npred)');
bayestopt_.restrict_columns = [ic; length(k2)+(1:nspred-npred)'];
aux = dr.transition_auxiliary_variables;
aux(:,2) = aux(:,2) + sum(k2 <= nstatic);
k = find(aux(:,2) > npred);
aux(k,2) = aux(k,2) + sum(k2 > nstatic+npred);
bayestopt_.restrict_aux = aux;


%% Initialization with unit-root variables
if ~isempty(options_.unit_root_vars)
  n_ur = size(options_.unit_root_vars,1);
  i_ur = zeros(n_ur,1);
  for i=1:n_ur
    i1 = strmatch(deblank(options_.unit_root_vars(i,:)),M_.endo_names(dr.order_var,:),'exact');
    if isempty(i1)
      error('Undeclared variable in unit_root_vars statement')
    end
    i_ur(i) = i1;
  end
  
  [junk,bayestopt_.var_list_stationary] = ...
      setdiff((1:M_.endo_nbr)',i_ur);
  [junk,bayestopt_.restrict_var_list_stationary] = ...
      setdiff(bayestopt_.restrict_var_list,i_ur);
  [junk,bayestopt_.restrict_var_list_nonstationary] = ...
      setdiff(bayestopt_.restrict_var_list,i_ur);
  if M_.maximum_lag > 1
    l1 = flipud([cumsum(M_.lead_lag_incidence(1:M_.maximum_lag-1,dr.order_var),1);ones(1,M_.endo_nbr)]);
    l2 = l1(:,restrict_var_list);
    il2 = find(l2' > 0);
    l2(il2) = (1:length(il2))';
    bayestopt_.restict_var_list_stationary = ...
	nonzeros(l2(:,restrict_var_list_stationary)); 
    bayestopt_.restict_var_list_nonstationary = ...
	nonzeros(l2(:,restrict_var_list_nonstationary)); 
  end
  options_.lik_init = 3;
end % if ~isempty(options_.unit_root_vars)

if isempty(options_.datafile)
  error('ESTIMATION: datafile option is missing')
end

%% If jscale isn't specified for an estimated parameter, use
%% global option options_.jscale, set to 0.2, by default
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

%% Read and demean data 
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);

options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1 & ~options_.logdata
  rawdata = log(rawdata);
end
if options_.prefilter == 1
  bayestopt_.mean_varobs = mean(rawdata,1);
  data = transpose(rawdata-ones(gend,1)*bayestopt_.mean_varobs);
else
  data = transpose(rawdata);
end

if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
	 ' transformation'])
end
if length(options_.mode_file) > 0 & options_.posterior_mode_estimation
  eval(['load ' options_.mode_file ';']');
end

%% compute sample moments if needed (bvar-dsge)
if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
  evalin('base',['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
                 'var_sample_moments(options_.first_obs,options_.first_obs+options_.nobs-1,options_.varlag,-1);'])
end

%% Compute the steadyn state if the _steadystate.m file is provided
if options_.steadystate_flag
  [oo_.steady_state,tchek] = feval([M_.fname '_steadystate'],[],[]);
end
initial_estimation_checks(xparam1,gend,data);

if options_.mode_compute == 0 & length(options_.mode_file) == 0
  return;
end


%% Estimation of the posterior mode or likelihood mode
if options_.mode_compute > 0 & options_.posterior_mode_estimation
  if isempty(strmatch('dsge_prior_weight',M_.param_names))
    fh=str2func('DsgeLikelihood');
  else
    fh=str2func('DsgeVarLikelihood');
  end
  if options_.mode_compute == 1
    optim_options = optimset('display','iter','LargeScale','off', ...
			     'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
    if isfield(options_,'optim_opt')
      eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if isempty(strmatch('dsge_prior_weight',M_.param_names))
      [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
          fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend,data);
    else
      [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
          fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend);
    end
  elseif options_.mode_compute == 2
    % asamin('set','maximum_cost_repeat',0);
    if isempty(strmatch('dsge_prior_weight',M_.param_names))
      [fval,xparam1,grad,hessian_asamin,exitflag] = ...
          asamin('minimize','DsgeLikelihood',xparam1,lb,ub,-ones(size(xparam1)),gend,data);   
    else
      [fval,xparam1,grad,hessian_asamin,exitflag] = ...
          asamin('minimize','DsgeVarLikelihood',xparam1,lb,ub,-ones(size(xparam1)),gend);   
    end       
  elseif options_.mode_compute == 3
    optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
    if isfield(options_,'optim_opt')
      eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if isempty(strmatch('dsge_prior_weight',M_.param_names))
      [xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend,data);
    else
      [xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend);
    end
  elseif options_.mode_compute == 4
    H0 = 1e-4*eye(nx);
    crit = 1e-7;
    nit = 1000;
    verbose = 2;
    if isempty(strmatch('dsge_prior_weight',M_.param_names))
      [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
          csminwel('DsgeLikelihood',xparam1,H0,[],crit,nit,gend,data);
      disp(sprintf('Objective function at mode: %f',fval))
      disp(sprintf('Objective function at mode: %f',DsgeLikelihood(xparam1,gend,data)))
    else
      [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
          csminwel('DsgeVarLikelihood',xparam1,H0,[],crit,nit,gend);
      disp(sprintf('Objective function at mode: %f',fval))
      disp(sprintf('Objective function at mode: %f',DsgeVarLikelihood(xparam1,gend)))
    end
  elseif options_.mode_compute == 5
    if isfield(options_,'hess')
      flag = options_.hess;
    else
      flag = 1;
    end
    if ~exist('igg'),  % by M. Ratto
      hh=[];
      gg=[];
      igg=[];
    end   % by M. Ratto
    if isfield(options_,'ftol')
      crit = options_.ftol;
    else
      crit = 1.e-7;
    end
    if isfield(options_,'nit')
      nit = options_.nit;
    else
      nit=1000;
    end
    if isempty(strmatch('dsge_prior_weight',M_.param_names))
      [xparam1,hh,gg,fval,invhess] = newrat('DsgeLikelihood',xparam1,hh,gg,igg,crit,nit,flag,gend,data);
    else
      [xparam1,hh,gg,fval,invhess] = newrat('DsgeVarLikelihood',xparam1,hh,gg,igg,crit,nit,flag,gend);
    end
    save([M_.fname '_mode'],'xparam1','hh','gg','fval','invhess');
    %eval(['save ' M_.fname '_mode xparam1 hh gg fval invhess;']);
  elseif options_.mode_compute == 6
    fval = DsgeLikelihood(xparam1,gend,data);
    OldMode = fval;
    if ~exist('MeanPar')
      MeanPar = xparam1;
    end
    if exist('hh')
      CovJump = inv(hh);
    else% The covariance matrix is initialized with the prior
        % covariance (a diagonal matrix) % Except for infinite variances ;-)
      stdev = bayestopt_.pstdev;
      indx = find(isinf(stdev));
      stdev(indx) = ones(length(indx),1)*0.1;
      indx = find(stdev>2);
      stdev(indx) = ones(length(indx),1)*0.1;      
      CovJump = diag(stdev).^2;
      CovJump = eye(length(stdev))*0.5;
    end
    OldPostVar = CovJump;
    Scale = options_.mh_jscale;
    for i=1:options_.Opt6Iter
      if i == 1
        if options_.Opt6Iter > 1
          flag = '';
        else
          flag = 'LastCall';
        end
        if isempty(strmatch('dsge_prior_weight',M_.param_names))
          [xparam1,PostVar,Scale,PostMean] = ...
              gmhmaxlik('DsgeLikelihood',xparam1,bounds,options_.Opt6Numb,Scale,flag,MeanPar,CovJump,gend,data);
          fval = DsgeLikelihood(xparam1,gend,data);
        else
          [xparam1,PostVar,Scale,PostMean] = ...
              gmhmaxlik('DsgeVarLikelihood',xparam1,bounds,options_.Opt6Numb,Scale,flag,MeanPar,CovJump,gend);
          fval = DsgeVarLikelihood(xparam1,gend);
        end
        options_.mh_jscale = Scale;
        mouvement = max(max(abs(PostVar-OldPostVar)));
        disp(['Change in the covariance matrix = ' num2str(mouvement) '.'])
        disp(['Mode improvement = ' num2str(abs(OldMode-fval))])
        OldMode = fval;
      else
        OldPostVar = PostVar;
        if i<options_.Opt6Iter
          flag = '';
        else
          flag = 'LastCall';
        end
        if isempty(strmatch('dsge_prior_weight',M_.param_names))
          [xparam1,PostVar,Scale,PostMean] = ...
              gmhmaxlik('DsgeLikelihood',xparam1,bounds,...
                        options_.Opt6Numb,Scale,flag,PostMean,PostVar,gend,data);
          fval = DsgeLikelihood(xparam1,gend,data);
        else
          [xparam1,PostVar,Scale,PostMean] = ...
              gmhmaxlik('DsgeVarLikelihood',xparam1,bounds,...
                        options_.Opt6Numb,Scale,flag,PostMean,PostVar,gend);
          fval = DsgeVarLikelihood(xparam1,gend);          
        end
        options_.mh_jscale = Scale;
        mouvement = max(max(abs(PostVar-OldPostVar)));
        fval = DsgeLikelihood(xparam1,gend,data);
        disp(['Change in the covariance matrix = ' num2str(mouvement) '.'])
        disp(['Mode improvement = ' num2str(abs(OldMode-fval))])
        OldMode = fval;
      end
      bayestopt_.jscale = ones(length(xparam1),1)*Scale;%??!
    end
    hh = inv(PostVar);
  end
  if options_.mode_compute ~= 5
    if options_.mode_compute ~= 6
      if isempty(strmatch('dsge_prior_weight',M_.param_names))
	hh = reshape(hessian('DsgeLikelihood',xparam1,gend,data),nx,nx);
      else
	hh = reshape(hessian('DsgeVarLikelihood',xparam1,gend),nx,nx);
      end
      save([M_.fname '_mode'],'xparam1','hh','fval');
      %eval(['save ' M_.fname '_mode xparam1 hh fval;']);
    else
      save([M_.fname '_mode'],'xparam1','hh','fval');
      %eval(['save ' M_.fname '_mode xparam1 hh fval;']);
    end
  end
  save([M_.fname '_mode'],'xparam1','hh');
  %eval(['save ' M_.fname '_mode xparam1 hh;']);
end

if options_.mode_check == 1 & options_.posterior_mode_estimation
  mode_check(xparam1,0,hh,gend,data,lb,ub);
end

if options_.posterior_mode_estimation
  hh = generalized_cholesky(hh);
  invhess = inv(hh);
  stdh = sqrt(diag(invhess));
else
  variances = bayestopt_.pstdev.^2;
  invhess = 0.001*diag(variances);
  invhess = 0.001*eye(length(variances));
end


if any(bayestopt_.pshape > 0) & options_.posterior_mode_estimation
  disp(' ')
  disp('RESULTS FROM POSTERIOR MAXIMIZATION')
  tstath = zeros(nx,1);
  for i = 1:nx
    tstath(i) = abs(xparam1(i))/stdh(i);
  end
  tit1 = sprintf('%10s %7s %8s %7s %6s %4s %6s\n',' ','prior mean', ...
		 'mode','s.d.','t-stat','prior','pstdev');
  if np
    ip = nvx+nvn+ncx+ncn+1;
    disp('parameters')
    disp(tit1)
    for i=1:np
      name = bayestopt_.name{ip};
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   name, ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mode.parameters.' name ' = xparam1(ip);']);
      eval(['oo_.posterior_std.parameters.' name ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if nvx
    ip = 1;
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:nvx
      k = estim_params_.var_exo(i,1);
      name = deblank(M_.exo_names(k,:));
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   name,bayestopt_.pmean(ip),xparam1(ip), ...
		   stdh(ip),tstath(ip),pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip))); 
      M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
      eval(['oo_.posterior_mode.shocks_std.' name ' = xparam1(ip);']);
      eval(['oo_.posterior_std.shocks_std.' name ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if nvn
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = nvx+1;
    for i=1:nvn
      name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   name,bayestopt_.pmean(ip), ...
		   xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mode.measurement_errors_std.' name ' = xparam1(ip);']);
      eval(['oo_.posterior_std.measurement_errors_std.' name ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if ncx
    disp('correlation of shocks')
    disp(tit1)
    ip = nvx+nvn+1;
    for i=1:ncx
      k1 = estim_params_.corrx(i,1);
      k2 = estim_params_.corrx(i,2);
      name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
      NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip),  ...
		   pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
      M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
      M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
      eval(['oo_.posterior_mode.shocks_corr.' NAME ' = xparam1(ip);']);
      eval(['oo_.posterior_std.shocks_corr.' NAME ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if ncn
    disp('correlation of measurement errors')
    disp(tit1)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
      NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
      eval(['oo_.posterior_mode.measurement_errors_corr.' NAME ' = xparam1(ip);']);
      eval(['oo_.posterior_std.measurement_errors_corr.' NAME ' = stdh(ip);']); 
      ip = ip+1;
    end
  end  
  %% Laplace approximation to the marginal log density:
  if isempty(strmatch('dsge_prior_weight',M_.param_names))
    md_Laplace = .5*size(xparam1,1)*log(2*pi) + .5*log(det(invhess)) ...
        - DsgeLikelihood(xparam1,gend,data);
  else
    md_Laplace = .5*size(xparam1,1)*log(2*pi) + .5*log(det(invhess)) ...
        - DsgeVarLikelihood(xparam1,gend);
  end
  oo_.MarginalDensity.LaplaceApproximation = md_Laplace;    
  disp(' ')
  disp(sprintf('Log data density [Laplace approximation] is %f.',md_Laplace))
  disp(' ')
elseif ~any(bayestopt_.pshape > 0) & options_.posterior_mode_estimation
  disp(' ')
  disp('RESULTS FROM MAXIMUM LIKELIHOOD')
  tstath = zeros(nx,1);
  for i = 1:nx
    tstath(i) = abs(xparam1(i))/stdh(i);
  end
  tit1 = sprintf('%10s %10s %7s %6s\n',' ','Estimate','s.d.','t-stat');
  if np
    ip = nvx+nvn+ncx+ncn+1;
    disp('parameters')
    disp(tit1)
    for i=1:np
      name = bayestopt_.name{ip};
      disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
		   name,xparam1(ip),stdh(ip),tstath(ip)));
      eval(['oo_.mle_mode.parameters.' name ' = xparam1(ip);']);
      eval(['oo_.mle_std.parameters.' name ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if nvx
    ip = 1;
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:nvx
      k = estim_params_.var_exo(i,1);
      name = deblank(M_.exo_names(k,:));
      disp(sprintf('%12s %8.4f %7.4f %7.4f',name,xparam1(ip),stdh(ip),tstath(ip)));
      M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
      eval(['oo_.mle_mode.shocks_std.' name ' = xparam1(ip);']);
      eval(['oo_.mle_std.shocks_std.' name ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if nvn
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = nvx+1;
    for i=1:nvn
      name = deblank(options_.varobs(estim_params_.var_endo(i,1),:));
      disp(sprintf('%12s %8.4f %7.4f %7.4f',name,xparam1(ip),stdh(ip),tstath(ip)))
      eval(['oo_.mle_mode.measurement_errors_std.' name ' = xparam1(ip);']);
      eval(['oo_.mle_std.measurement_errors_std.' name ' = stdh(ip);']);      
      ip = ip+1;
    end
  end
  if ncx
    disp('correlation of shocks')
    disp(tit1)
    ip = nvx+nvn+1;
    for i=1:ncx
      k1 = estim_params_.corrx(i,1);
      k2 = estim_params_.corrx(i,2);
      name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
      NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
      disp(sprintf('%12s %8.4f %7.4f %7.4f', name,xparam1(ip),stdh(ip),tstath(ip)));
      M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
      M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
      eval(['oo_.mle_mode.shocks_corr.' NAME ' = xparam1(ip);']);
      eval(['oo_.mle_std.shocks_corr.' NAME ' = stdh(ip);']);      
      ip = ip+1;
    end
  end
  if ncn
    disp('correlation of measurement errors')
    disp(tit1)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
      NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
      disp(sprintf('%12s %8.4f %7.4f %7.4f',name,xparam1(ip),stdh(ip),tstath(ip)));
      eval(['oo_.mle_mode.measurement_error_corr.' NAME ' = xparam1(ip);']);
      eval(['oo_.mle_std.measurement_error_corr.' NAME ' = stdh(ip);']);
      ip = ip+1;
    end
  end
end


OutputDirectoryName = CheckPath('Output');

if any(bayestopt_.pshape > 0) & options_.TeX %% Bayesian estimation (posterior mode) Latex output
  if np
    filename = [OutputDirectoryName '\' M_.fname '_Posterior_Mode_1.TeX'];
    fidTeX = fopen(filename,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (parameters)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'{\\tiny \n')
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+ncx+ncn+1;
    for i=1:np
      fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      M_.param_names_tex(estim_params_.param_vals(i,1),:),...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
	      bayestopt_.pmean(ip),...
	      estim_params_.param_vals(i,6),...
	      xparam1(ip),...
	      stdh(ip));
      ip = ip + 1;    
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (parameters)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:1}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'} \n')
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if nvx
    TeXfile = [OutputDirectoryName '\' M_.fname '_Posterior_Mode_2.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of structural shocks)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'{\\tiny \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = 1;
    for i=1:nvx
      k = estim_params_.var_exo(i,1);
      fprintf(fidTeX,[ '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
	      deblank(M_.exo_names_tex(k,:)),...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
	      bayestopt_.pmean(ip),...
	      estim_params_.var_exo(i,7),...
	      xparam1(ip), ...
	      stdh(ip)); 
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of structural shocks)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:2}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'} \n')
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if nvn
    TeXfile = [OutputDirectoryName '\' M_.fname '_Posterior_Mode_3.TeX'];
    fidTeX  = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of measurement errors)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+1;
    for i=1:nvn
      idx = strmatch(options_.varobs(estim_params_.var_endo(i,1),:),M_.endo_names);
      fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      deblank(M_.endo_names_tex(idx,:)), ...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...        
	      bayestopt_.pmean(ip), ...
	      estim_params_.var_endo(i,7),...        
	      xparam1(ip),...
	      stdh(ip)); 
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of measurement errors)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:3}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if ncx
    TeXfile = [OutputDirectoryName '\' M_.fname '_Posterior_Mode_4.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of structural shocks)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+1;
    for i=1:ncx
      k1 = estim_params_.corrx(i,1);
      k2 = estim_params_.corrx(i,2);
      fprintf(fidTeX,[ '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
	      [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))], ...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
	      bayestopt_.pmean(ip), ...
	      estim_params_.corrx(i,8), ...
	      xparam1(ip), ...
	      stdh(ip));
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of structural shocks)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:4}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if ncn
    TeXfile = [OutputDirectoryName '\' M_.fname '_Posterior_Mode_5.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of measurement errors)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))], ...
	      pnames(bayestopt_.pshape(ip)+1,:), ...
	      bayestopt_.pmean(ip), ...
	      estim_params_.corrn(i,8), ...
	      xparam1(ip), ...
	      stdh(ip));
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of measurement errors)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:5}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
end

pindx = estim_params_.param_vals(:,1);
save([M_.fname '_params'],'pindx');

if (any(bayestopt_.pshape  >0 ) & options_.mh_replic) | ...
      (any(bayestopt_.pshape >0 ) & options_.load_mh_file)  %% not ML estimation
  bounds = prior_bounds(bayestopt_);
  bayestopt_.lb = bounds(:,1);
  bayestopt_.ub = bounds(:,2);
  if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
    find(xparam1 < bounds(:,1))
    find(xparam1 > bounds(:,2))
    error('Mode values are outside prior bounds. Reduce prior_trunc.')
  end
  if options_.mh_replic
    if ~options_.load_mh_file
      if isempty(strmatch('dsge_prior_weight',M_.param_names))
        metropolis('DsgeLikelihood',xparam1,invhess,bounds,gend,data);
      else
        metropolis('DsgeVarLikelihood',xparam1,invhess,bounds,gend);
      end
    else
      if options_.use_mh_covariance_matrix
        invhess = compute_mh_covariance_matrix();
      end
      if isempty(strmatch('dsge_prior_weight',M_.param_names))
        metropolis('DsgeLikelihood',xparam1,invhess,bounds,gend,data);
      else
        metropolis('DsgeVarLikelihood',xparam1,invhess,bounds,gend);
      end
    end
  end
  if ~options_.nodiagnostic & options_.mh_replic > 1000 & options_.mh_nblck > 1
    McMCDiagnostics;
  end
  %% Here i discard first half of the draws:
  CutSample;
  %% Estimation of the marginal density from the Mh draws:
  marginal = marginal_density;
  %% 
  GetPosteriorParametersStatistics;
  PlotPosteriorDistributions;
  metropolis_draw(1);
  if options_.bayesian_irf
    PosteriorIRF('posterior');
  end
  return
end

if ~((any(bayestopt_.pshape > 0) & options_.mh_replic) | (any(bayestopt_.pshape ...
						  > 0) & options_.load_mh_file)) | ~options_.smoother  
    %% ML estimation, or posterior mode without metropolis-hastings or metropolis without bayesian smooth variables
  options_.lik_algo = 2;
  [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(xparam1,gend,data);
  for i=1:M_.endo_nbr
    eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(dr.order_var(i),:)) ' = atT(i,:)'';']);
    eval(['oo_.FilteredVariables.' deblank(M_.endo_names(dr.order_var(i),:)) ' = filtered_state_vector(i,:)'';']);
  end
  [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
  if options_.TeX
    fidTeX = fopen([M_.fname '_SmoothedShocks.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
  end    
  if nbplt == 1
    hh = figure('Name','Smoothed shocks');
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:M_.exo_nbr
      subplot(nr,nc,i);
      plot(1:gend,innov(i,:),'-k','linewidth',1)
      hold on
      plot([1 gend],[0 0],'-r','linewidth',.5)
      hold off
      xlim([1 gend])
      name    = deblank(M_.exo_names(i,:));
      NAMES   = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
	set(gca,'XTick',options_.XTick)
	set(gca,'XTickLabel',options_.XTickLabel)
      end
      if options_.TeX
	texname = M_.exo_names_tex(i,1);
	TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none')
      eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
    end
    eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(1)]);
    eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(1)]);
    saveas(hh,[M_.fname '_SmoothedShocks' int2str(1) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:M_.exo_nbr
	fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(1));
      fprintf(fidTeX,'\\caption{Smoothed shocks.}');
      fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(1));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end
  else
    for plt = 1:nbplt-1
      hh = figure('Name','Smoothed shocks');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:nstar
	k = (plt-1)*nstar+i;
	subplot(nr,nc,i);
	plot([1 gend],[0 0],'-r','linewidth',.5)
	hold on
	plot(1:gend,innov(k,:),'-k','linewidth',1)
	hold off
	name = deblank(M_.exo_names(k,:));
	NAMES = strvcat(NAMES,name);
	if ~isempty(options_.XTick)
	  set(gca,'XTick',options_.XTick)
	  set(gca,'XTickLabel',options_.XTickLabel)
	end
	xlim([1 gend])
	if options_.TeX
	  texname = M_.exo_names_tex(k,:);
	  TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
	end    
	title(name,'Interpreter','none')
	eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
      end
      eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(plt)]);
      eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(plt)]);
      saveas(hh,[M_.fname '_SmoothedShocks' int2str(plt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
	fprintf(fidTeX,'\\begin{figure}[H]\n');
	for jj = 1:nstar
	  fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
	end    
	fprintf(fidTeX,'\\centering \n');
	fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(plt));
	fprintf(fidTeX,'\\caption{Smoothed shocks.}');
	fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
	fprintf(fidTeX,'\\end{figure}\n');
	fprintf(fidTeX,'\n');
      end    
    end
    hh = figure('Name','Smoothed shocks');
    set(0,'CurrentFigure',hh)
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:M_.exo_nbr-(nbplt-1)*nstar
      k = (nbplt-1)*nstar+i;
      if lr ~= 0
        subplot(lr,lc,i);
      else
        subplot(nr,nc,i);
      end    
      plot([1 gend],[0 0],'-r','linewidth',0.5)
      hold on
      plot(1:gend,innov(k,:),'-k','linewidth',1)
      hold off
      name     = deblank(M_.exo_names(k,:));
      NAMES    = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      xlim([1 gend])
      if options_.TeX
        texname  = M_.exo_names_tex(k,:);
        TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none')
      eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
    end
    eval(['print -depsc2 ' M_.fname '_SmoothedShocks' int2str(nbplt)]);
    eval(['print -dpdf ' M_.fname '_SmoothedShocks' int2str(nbplt)]);
    saveas(hh,[M_.fname '_SmoothedShocks' int2str(nbplt) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:size(NAMES,1);
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(nbplt));
      fprintf(fidTeX,'\\caption{Smoothed shocks.}');
      fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(nbplt));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
  end
  %%
  %%	Smooth observational errors...
  %%
  yf = zeros(gend,n_varobs);
  if options_.prefilter == 1
    yf = atT(bayestopt_.mf,:)+repmat(transpose(bayestopt_.mean_varobs),1,gend);
  elseif options_.loglinear == 1
    yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
         trend_coeff*[1:gend];
  else
    yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
         trend_coeff*[1:gend];
  end
  if nvn
    number_of_plots_to_draw = 0;
    index = [];
    for i=1:n_varobs
      if max(abs(measurement_error(10:end))) > 0.000000001
        number_of_plots_to_draw = number_of_plots_to_draw + 1;
        index = cat(1,index,i);
      end
      eval(['oo_.SmoothedMeasurementErrors.' deblank(options_.varobs(i,:)) ...
            ' = measurement_error(i,:)'';']);
    end
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
    if options_.TeX
      fidTeX = fopen([M_.fname '_SmoothedObservationErrors.TeX'],'w');
      fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
      fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      fprintf(fidTeX,' \n');
    end
    if nbplt == 1
      hh = figure('Name','Smoothed observation errors');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:number_of_plots_to_draw
        subplot(nr,nc,i);
        plot(1:gend,measurement_error(index(i),:),'-k','linewidth',1)
        hold on
        plot([1 gend],[0 0],'-r','linewidth',.5)
        hold off
        name    = deblank(options_.varobs(index(i),:));
        NAMES   = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        if options_.TeX
          idx = strmatch(options_.varobs(indx(i),:),M_.endo_names,'exact');
          texname = M_.endo_names_tex(idx,:);
          TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end
        title(name,'Interpreter','none')
      end
      eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(1)]);
      eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(1)]);
      saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(1) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:number_of_plots_to_draw
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(1));
        fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
        fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s',int2str(1));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
      end
    else
      for plt = 1:nbplt-1
        hh = figure('Name','Smoothed observation errors');
        set(0,'CurrentFigure',hh)
        NAMES = [];
        if options_.TeX, TeXNAMES = [], end
        for i=1:nstar
          k = (plt-1)*nstar+i;
          subplot(nr,nc,i);
          plot([1 gend],[0 0],'-r','linewidth',.5)
          hold on
          plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
          hold off
          name = deblank(options_.varobs(index(k),:));
          NAMES = strvcat(NAMES,name);
          if ~isempty(options_.XTick)
            set(gca,'XTick',options_.XTick)
            set(gca,'XTickLabel',options_.XTickLabel)
          end
          if options_.TeX
            idx = strmatch(options_.varobs(k),M_.endo_names,'exact');
            texname = M_.endo_names_tex(idx,:);
            TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
          end    
          title(name,'Interpreter','none')
        end
        eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(plt)]);
        eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(plt)]);
        saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(plt) '.fig']);
        if options_.nograph, close(hh), end
        if options_.TeX
          fprintf(fidTeX,'\\begin{figure}[H]\n');
          for jj = 1:nstar
            fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
          end    
          fprintf(fidTeX,'\\centering \n');
          fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(plt));
          fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
          fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
          fprintf(fidTeX,'\\end{figure}\n');
          fprintf(fidTeX,'\n');
        end    
      end
      hh = figure('Name','Smoothed observation errors');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:number_of_plots_to_draw-(nbplt-1)*nstar
        k = (nbplt-1)*nstar+i;
        if lr ~= 0
          subplot(lr,lc,i);
        else
          subplot(nr,nc,i);
        end    
        plot([1 gend],[0 0],'-r','linewidth',0.5)
        hold on
        plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
        hold off
        name     = deblank(options_.varobs(index(k),:));
        NAMES    = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        if options_.TeX
          idx = strmatch(options_.varobs(index(k)),M_.endo_names,'exact');
          texname = M_.endo_names_tex(idx,:);
          TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end
        title(name,'Interpreter','none');
      end
      eval(['print -depsc2 ' M_.fname '_SmoothedObservationErrors' int2str(nbplt)]);
      eval(['print -dpdf ' M_.fname '_SmoothedObservationErrors' int2str(nbplt)]);
      saveas(hh,[M_.fname '_SmoothedObservationErrors' int2str(nbplt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:size(NAMES,1);
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservedErrors%s}\n',M_.fname,int2str(nbplt));
        fprintf(fidTeX,'\\caption{Smoothed observed errors.}');
        fprintf(fidTeX,'\\label{Fig:SmoothedObservedErrors:%s}\n',int2str(nbplt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
      end    
    end
  end	
  %%
  %%	Historical and smoothed variabes
  %%
  [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
  if options_.TeX
    fidTeX = fopen([M_.fname '_HistoricalAndSmoothedVariables.TeX'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
  end    
  if nbplt == 1
    hh = figure('Name','Historical and smoothed variables');
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:n_varobs
      subplot(nr,nc,i);
      plot(1:gend,yf(i,:),'-r','linewidth',1)
      hold on
      plot(1:gend,rawdata(:,i),'-k','linewidth',1)
      hold off
      name    = deblank(options_.varobs(i,:));
      NAMES   = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      xlim([1 gend])
      if options_.TeX
        idx = strmatch(options_.varobs(i),M_.endo_names,'exact');
        texname = M_.endo_names_tex(idx,:);
        TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none')
    end
    eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1)]);
    eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(1)]);
    saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(1) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:n_varobs
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(1));
      fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
      fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(1));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
  else
    for plt = 1:nbplt-1
      hh = figure('Name','Historical and smoothed variables');
      set(0,'CurrentFigure',hh)
      NAMES = [];
      if options_.TeX, TeXNAMES = [], end
      for i=1:nstar
        k = (plt-1)*nstar+i;
        subplot(nr,nc,i);
        plot(1:gend,yf(k,:),'-r','linewidth',1)
        hold on
        plot(1:gend,rawdata(:,k),'-k','linewidth',1)
        hold off
        name = deblank(options_.varobs(k,:));
        NAMES = strvcat(NAMES,name);
        if ~isempty(options_.XTick)
          set(gca,'XTick',options_.XTick)
          set(gca,'XTickLabel',options_.XTickLabel)
        end
        xlim([1 gend])
        if options_.TeX
          idx = strmatch(options_.varobs(k),M_.endo_names,'exact');
          texname = M_.endo_names_tex(idx,:);
          TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
        end    
        title(name,'Interpreter','none')
      end
      eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)]);
      eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)]);
      saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(plt) '.fig']);
      if options_.nograph, close(hh), end
      if options_.TeX
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        for jj = 1:nstar
          fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
        end    
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(plt));
        fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
        fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
      end    
    end
    hh = figure('Name','Historical and smoothed variables');
    set(0,'CurrentFigure',hh)
    NAMES = [];
    if options_.TeX, TeXNAMES = [], end
    for i=1:n_varobs-(nbplt-1)*nstar
      k = (nbplt-1)*nstar+i;
      if lr ~= 0
        subplot(lr,lc,i);
      else
        subplot(nr,nc,i);
      end    
      plot(1:gend,yf(k,:),'-r','linewidth',1)
      hold on
      plot(1:gend,rawdata(:,k),'-k','linewidth',1)
      hold off
      name = deblank(options_.varobs(k,:));
      NAMES    = strvcat(NAMES,name);
      if ~isempty(options_.XTick)
        set(gca,'XTick',options_.XTick)
        set(gca,'XTickLabel',options_.XTickLabel)
      end
      xlim([1 gend])
      if options_.TeX
        idx = strmatch(options_.varobs(i),M_.endo_names,'exact');
        texname = M_.endo_names_tex(idx,:);
        TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      end
      title(name,'Interpreter','none');
    end
    eval(['print -depsc2 ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
    eval(['print -dpdf ' M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
    saveas(hh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(nbplt) '.fig']);
    if options_.nograph, close(hh), end
    if options_.TeX
      fprintf(fidTeX,'\\begin{figure}[H]\n');
      for jj = 1:size(NAMES,1);
        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      end    
      fprintf(fidTeX,'\\centering \n');
      fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(nbplt));
      fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
      fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(nbplt));
      fprintf(fidTeX,'\\end{figure}\n');
      fprintf(fidTeX,'\n');
      fprintf(fidTeX,'%% End of TeX file.\n');
      fclose(fidTeX);
    end    
  end
end 

if options_.forecast > 0 & options_.mh_replic == 0 & ~options_.load_mh_file 
  forecast(var_list);
end

pindx = estim_params_.param_vals(:,1);
save([M_.fname '_pindx','pindx']);