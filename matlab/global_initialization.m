function global_initialization()
% initializes global variables and options for DYNARE
  global oo_ M_ options_ ct_ endval_ rplottype_
  
  ct_=0;
  endval_=0;
  
  options_.rplottype = 0;
  options_.smpl = 0;
  options_.dynatol = 0.00001;
  options_.maxit_ = 10;
  options_.slowc = 1;
  options_.timing = 0;
  options_.gstep = 1e-2;
  options_.debug = 0;
  options_.initval_file = 0;
  options_.Schur_vec_tol = 1e-8; % used to find nonstationary variables
                                 % in Schur decomposition of the
                                 % transition matrix
  options_.solve_tolf = eps^(1/3);
  options_.solve_tolx = eps^(2/3);
  options_.solve_maxit = 500;


  % steady state file
  if exist([M_.fname '_steadystate'])
    options_.steadystate_flag = 1;
  else
    options_.steadystate_flag = 0;
  end
  options_.steadystate_partial = [];

  % subset of the estimated deep parameters
  options_.ParamSubSet = 'None';

  % bvar-dsge
  options_.bvar_dsge = 0;
  options_.varlag = 4;

  % Optimization algorithm [6] gmhmaxlik
  options_.Opt6Iter = 3;
  options_.Opt6Numb = 100000;

  % Graphics  
  options_.graphics.nrows = 3;
  options_.graphics.ncols = 3;
  options_.graphics.line_types = {'b-'};
  options_.graphics.line_width = 1;
  options_.nograph = 0;
  options_.XTick = [];
  options_.XTickLabel = [];

  % IRFs & other stoch_simul output
  options_.irf = 40;
  options_.relative_irf = 0;
  options_.ar = 5;
  options_.simul_seed = [];
  options_.hp_filter = 0;
  options_.hp_ngrid = 512;
  options_.nomoments = 0;
  options_.nocorr = 0;
  options_.periods = 0;
  options_.noprint = 0;
  options_.simul = 0;
  options_.SpectralDensity = 0;

  % TeX output
  options_.TeX = 0;

  % Exel
  options_.xls_sheet = '';
  options_.xls_range = '';

  % Prior draws
  options_.forecast = 0;
  options_.replic = 1;

  % Model
  options_.linear = 0;

  % Solution
  options_.order = 2;
  options_.dr_algo = 0;
  options_.solve_algo = 2;
  options_.linear = 0;
  options_.replic = 50;
  options_.drop = 100;
  options_.simul_algo = 0;

  % Ramsey policy
  options_.planner_discount = 1.0;
  options_.ramsey_policy = 0;

  % estimation
  options_.load_mh_file = 0;
  options_.first_obs = 1;
  options_.prefilter = 0;
  options_.presample = 0;
  options_.lik_algo = 1;
  options_.lik_init = 1;
  options_.mh_replic = 20000;
  options_.mh_drop = 0.5;
  options_.mh_jscale = 0.2;
  options_.mh_init_scale = 2*options_.mh_jscale;
  options_.mode_file = '';
  options_.mode_compute = 4;
  options_.mode_check = 0;
  options_.prior_trunc = 1e-10;
  options_.mh_conf_sig = 0.90;
  options_.mh_mode = 1;
  options_.mh_nblck = 2;
  options_.load_mh_file = 0;
  options_.mh_recover = 0;
  options_.nodiagnostic = 0;
  options_.loglinear = 0;
  options_.unit_root_vars = [];
  options_.bayesian_irf = 0;
  options_.bayesian_th_moments = 0;
  options_.smoother = 0;
  options_.moments_varendo = 0;
  options_.filtered_vars = 0;
  options_.kalman_algo = 1;
  options_.kalman_tol = 1e-12;
  options_.posterior_mode_estimation = 1;
  options_.MaxNumberOfBytes = 1e6;
  options_.filter_step_ahead = [];
  options_.diffuse_d = [];
  options_.logdata = 0;
  options_.use_mh_covariance_matrix = 0;
  options_.noconstant = 0;
  options_.markowitz = 0.5;

  % Misc
  options_.conf_sig = 0.9;
  oo_.exo_simul = [];
  oo_.endo_simul = [];
  oo_.dr = [];
  oo_.exo_steady_state = [];
  oo_.exo_det_steady_state = [];
  oo_.exo_det_simul = [];

  % Variance matrix for measurement errors
  M_.H = 0;

  % BVAR
  M_.bvar = [];
