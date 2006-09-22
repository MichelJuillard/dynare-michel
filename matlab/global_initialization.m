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
  options_.debug = 0
  options_.initval_file = 0;
  options_.Schur_vec_tol = 1e-8; % used to find nonstationary variables
                                 % in Schur decomposition of the
                                 % transition matrix
  options_.solve_tolf = eps^(2/3);
  options_.solve_tolx = 3.7e-11;
  options_.solve_maxit = 500;
				 
  if exist([M_.fname '_steadystate'])
    options_.steadystate_flag = 1;
  else
    options_.steadystate_flag = 0;
  end
  options_.steadystate_partial = [];
  
  options_.ParamSubSet = 'None';
  
  % bvar-dsge
  options_.varlag = 4;
  
  % Optimization algorithm [6] gmhmaxlik
  options_.Opt6Iter = 3;
  options_.Opt6Numb = 100000;
  
  
  oo_.exo_simul = [];
  oo_.endo_simul = [];
  oo_.dr = struct([]);
  oo_.exo_det_steady_state = [];
  oo_.exo_det_simul = [];