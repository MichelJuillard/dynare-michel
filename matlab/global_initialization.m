function global_initialization()
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

  oo_.exo_simul = [];
  oo_.endo_simul = [];
  oo_.dr = struct([]);
  oo_.exo_det_steady_state = [];
  oo_.exo_det_simul = [];