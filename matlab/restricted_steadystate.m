function [sR,sG] = restricted_steadystate(y,x,indx)
  % stephane.adjemian@gmail.com
  global options_ M_ oo_
  
  inde  = options_.steadystate_partial.sseqn;
  
  ss = oo_.steady_state;
  
  ss(indx) = y;
 
  eval(['[R,G] = ' M_.fname '_static(ss, x, M_.params);']);

  sR = R(inde);
  sG = G(inde,indx);