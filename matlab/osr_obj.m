function [loss,vx,info]=osr_obj(x,i_params,weights);
  % objective function for optimal simple rules (OSR)
  global M_ oo_ optimal_Q_ it_
%  global ys_ Sigma_e_ endo_nbr exo_nbr optimal_Q_ it_ ykmin_ options_
  
  vx = [];
  % set parameters of the policiy rule
  M_.params(i_params) = x;
  
  % don't change below until the part where the loss function is computed
  it_ = M_.maximum_lag+1;
  [dr,info] = resol(oo_.steady_state,0);
  
  switch info(1)
   case 1
    loss = 1e8;
    return
   case 2
    loss = 1e8*min(1e3,info(2));
    return
   case 3
    loss = 1e8*min(1e3,info(2));
    return
   case 4
    loss = 1e8*min(1e3,info(2));
    return
   case 5
    loss = 1e8;
    return
   case 20
    loss = 1e8*min(1e3,info(2));
    return
   otherwise
  end
  
  [A,B] = kalman_transition_matrix(dr);
  [vx,ns_var] = lyapunov_symm(A,B*M_.Sigma_e*B');
  endo_nbr = M_.endo_nbr;
  i_var = (1:endo_nbr)';
  i_var(ns_var) = zeros(length(ns_var),1);
  i_var = nonzeros(i_var);
  vx = vx(i_var,i_var);
  weights = weights(dr.order_var,dr.order_var);
  weights = sparse(weights(i_var,i_var));
  
  loss = weights(:)'*vx(:);
  














