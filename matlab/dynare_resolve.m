function [A,B,ys] = dynare_resolve()
global dr1_test_ oo_

dr = resol1(oo_.steady_state,0,1,1);
oo_.dr = dr;

if dr1_test_(1) > 0
	A = [];
    B = [];
    ys = [];
    return
end
  
[A,B] = kalman_transition_matrix(dr);
ys = dr.ys;

% SA 01-18-2005 	Added line five (oo_.dr = dr;)
