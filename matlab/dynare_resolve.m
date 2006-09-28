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
    nstatic = oo_.dr.nstatic;
    npred = oo_.dr.npred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:npred) endo_nbr+(1:size(oo_.dr.ghx,2)-npred) ]';
    aux = oo_.dr.transition_auxiliary_variables;
    k = find(aux(:,2) > npred);
    aux(:,2) = aux(:,2) + nstatic;
    aux(k,2) = aux(k,2) + oo_.dr.nfwrd;
  end
  
  [A,B] = kalman_transition_matrix(oo_.dr,iv,ic,aux);
  ys = oo_.dr.ys;