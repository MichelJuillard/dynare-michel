function [vx1,i_ns] = get_variance_of_endogenous_variables(dr,i_var)

  global M_ options_
  
  Sigma_e = M_.Sigma_e;
  
  nstatic = dr.nstatic;
  npred = dr.npred;
  ghx = dr.ghx(i_var,:);
  ghu = dr.ghu(i_var,:);
  nc = size(ghx,2);
  n = length(i_var);
  
  [A,B] = kalman_transition_matrix(dr,nstatic+(1:npred),1:nc,dr.transition_auxiliary_variables);
  
  [vx,u] = lyapunov_symm(A,B*Sigma_e*B');
  
  if size(u,2) > 0
    i_stat = find(any(abs(ghx*u) < options_.Schur_vec_tol,2));
  
    ghx = ghx(i_stat,:);
    ghu = ghu(i_stat,:);
  else
    i_stat = (1:n)';
  end
  
  vx1 = Inf*ones(n,n);
  vx1(i_stat,i_stat) = ghx*vx*ghx'+ghu*Sigma_e*ghu';
  