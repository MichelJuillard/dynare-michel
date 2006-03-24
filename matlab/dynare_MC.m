function dynare_MC(var_list_)

global M_ options_ oo_ estim_params_ 
global bayestopt_

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
options_ = set_default_option(options_,'irf',0);
options_ = set_default_option(options_,'relative_irf',0);
options_ = set_default_option(options_,'order',1);
options_ = set_default_option(options_,'ar',5);
options_ = set_default_option(options_,'dr_algo',0);
options_ = set_default_option(options_,'linear',0);
options_ = set_default_option(options_,'drop',0);
options_ = set_default_option(options_,'replic',1);
options_ = set_default_option(options_,'hp_filter',0);
options_ = set_default_option(options_,'forecast',0);
options_ = set_default_option(options_,'smoother',0);
options_ = set_default_option(options_,'moments_varendo',0);
options_ = set_default_option(options_,'filtered_vars',0);
options_ = set_default_option(options_,'kalman_algo',1);
options_ = set_default_option(options_,'kalman_tol',10^(-12));
options_ = set_default_option(options_,'posterior_mode_estimation',1);
options_ = set_default_option(options_,'MaxNumberOfBytes',1e6);
options_ = set_default_option(options_,'xls_sheet','');
options_ = set_default_option(options_,'xls_range','');

%% Add something to the parser ++>
M_.dname = M_.fname; % The user should be able to choose another name
                     % for the directory...


pnames 		= ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
n_varobs 	= size(options_.varobs,1);

[xparam1,estim_params_,bayestopt_,lb,ub] = set_prior(estim_params_);

options_.mh_replic = 0;

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

%% Static solver
if exist([M_.fname '_steadystate'])
  bayestopt_.static_solve = [M_.fname '_steadystate'];
else
  bayestopt_.static_solve = 'dynare_solve';
end

dr = set_state_space([]);

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
  if M_.maximum_lag > 1
    l1 = flipud([cumsum(M_.lead_lag_incidence(1:M_.maximum_lag-1,dr.order_var),1);ones(1,M_.endo_nbr)]);
    n1 = nnz(l1);
    bayestopt_.Pinf = zeros(n1,n1);
    l2 = find(l1');
    l3 = zeros(M_.endo_nbr,M_.maximum_lag);
    l3(i_ur,:) = l1(:,i_ur)';
    l3 = l3(:);
    i_ur1 = find(l3(l2));
    i_stable = ones(M_.endo_nbr,1);
    i_stable(i_ur) = zeros(n_ur,1);
    i_stable = find(i_stable);
    bayestopt_.Pinf(i_ur1,i_ur1) = diag(ones(1,length(i_ur1)));
    bayestopt_.i_var_stable = i_stable;
    l3 = zeros(M_.endo_nbr,M_.maximum_lag);
    l3(i_stable,:) = l1(:,i_stable)';
    l3 = l3(:);
    bayestopt_.i_T_var_stable = find(l3(l2));
  else
    n1 = M_.endo_nbr;
    bayestopt_.Pinf = zeros(n1,n1);
    bayestopt_.Pinf(i_ur,i_ur) = diag(ones(1,length(i_ur)));
    l1 = ones(M_.endo_nbr,1);
    l1(i_ur,:) = zeros(length(i_ur),1);
    bayestopt_.i_T_var_stable = find(l1);
  end
  options_.lik_init = 3;
end % if ~isempty(options_.unit_root_vars)

if isempty(options_.datafile)
  error('ESTIMATION: datafile option is missing')
end

if isempty(options_.varobs)
  error('ESTIMATION: VAROBS is missing')
end


%% If jscale isn't specified for an estimated parameter, use
%% global option options_.jscale, set to 0.2, by default
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

%% Read and demean data 
rawdata = read_variables(options_.datafile,options_.varobs,[],options_.xls_sheet,options_.xls_range);

k = [];
k1 = [];
for i=1:n_varobs
  k = [k strmatch(deblank(options_.varobs(i,:)),M_.endo_names(dr.order_var,:), ...
		  'exact')];
  k1 = [k1 strmatch(deblank(options_.varobs(i,:)),M_.endo_names, 'exact')];
end

bayestopt_.mf 	= k;
bayestopt_.mfys = k1;
options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1
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

nvx  = estim_params_.nvx;
nvn  = estim_params_.nvn;
ncx  = estim_params_.ncx;
ncn  = estim_params_.ncn;
np   = estim_params_.np ;
npar = nvx+nvn+ncx+ncn+np;
offset = npar-np;
fname_=M_.fname;

options_ = set_default_option(options_,'opt_gsa',1);
options_gsa_ = options_.opt_gsa;

if options_gsa_.pprior,
  load([fname_,'_prior'])
  namfile=[fname_,'_prior'];
else
  load([fname_,'_mc'])
  namfile=[fname_,'_mc'];
end
load(options_.mode_file)
%%
%%
%%
x=[lpmat0(istable,:) lpmat(istable,:)];
clear lpmat lpmat0 istable iunstable egg yys T
B = size(x,1);
if M_.maximum_lag > 1
    l1 = flipud([cumsum(M_.lead_lag_incidence(1:M_.maximum_lag-1,dr.order_var),1);ones(1,M_.endo_nbr)]);
    n1 = nnz(l1);
else
    n1 = M_.endo_nbr;
end
stock_smooth = zeros(gend,n1,B);
stock_filter = zeros(gend+1,n1,B);
stock_ys = zeros(M_.endo_nbr,B);
%%
h = waitbar(0,'MC smoother ...');
for b=1:B
  deep = x(b,:);
  %deep(1:offset) = xparam1(1:offset);
  logpo2(b) = DsgeLikelihood(deep',gend,data);
  [atT,innov,measurement_error,filtered_state_vector,ys,trend_coeff] = DsgeSmoother(deep,gend,data);
  stock_smooth(:,:,b)=atT';
  stock_filter(:,:,b)=filtered_state_vector';
  stock_ys(:,b)=ys;
  waitbar(b/B,h,['MC smoother ...',num2str(b),'/',num2str(B)]);
end
close(h)

save(namfile,'x','logpo2','stock_smooth','stock_filter','stock_ys','-append')
