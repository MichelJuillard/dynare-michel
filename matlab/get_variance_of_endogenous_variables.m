function [vx,i_ns] = get_variance_of_endogenous_variables(dr,i_var)

  global M_ options_
  
  Sigma_e = M_.Sigma_e;
  
  nstatic = dr.nstatic;
  npred = dr.npred;
  ghx = dr.ghx(i_var,:);
  ghu = dr.ghu(i_var,:);
  nc = size(ghx,2);
  
  [A,B] = kalman_transition_matrix(dr,nstatic+(1:npred),1:nc,dr.transition_auxiliary_variables);
  
  [vx,u] = lyapunov_symm(A,B*Sigma_e*B');
  
  i_ns = find(any(abs(ghx*u) > options_.Schur_vec_tol,2));
  
  ghx = ghx(i_ns,:);
  ghu = ghu(i_ns,:);
  
  n = length(i_var);
  vx = Inf*ones(n,n);
  vx(i_ns,i_ns) = ghx*vx*ghx'+ghu*Sigma_e*ghu';
  