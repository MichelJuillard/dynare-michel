function [A,B,ys,info] = dynare_resolve(iv,ic,aux)
  global oo_ M_
  
  [oo_.dr,info] = resol(oo_.steady_state,0);
  
  if info(1) > 0
    A = [];
    B = [];
    ys = [];
    return
  end
  
  if nargin == 0
    endo_nbr = M_.endo_nbr;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:npred) endo_nbr+(1:size(dr.ghx,2)-npred) ]';
    aux = oo_.dr.transtion_auxiliary_variables(:,2);
    k = find(aux > oo_.dr.npred);
    aux = aux + nstatic;
    aux(k) = aux(k) + oo_.dr.nfrwd;
  end
  
  [A,B] = kalman_transition_matrix(oo_.dr,iv,ic,aux);
  ys = oo_.dr.ys;