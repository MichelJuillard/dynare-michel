function [A,B,ys,info] = dynare_resolve()
  global oo_
  
  [oo_.dr,info] = resol(oo_.steady_state,0);
  
  if info(1) > 0
    A = [];
    B = [];
    ys = [];
    return
  end
  
  [A,B] = kalman_transition_matrix(oo_.dr);
  ys = oo_.dr.ys;