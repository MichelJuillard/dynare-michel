function tau = thet2tau(params, indx, flagmoments,mf)
global M_ oo_ options_

if nargin==1,
    indx = [1:M_.param_nbr];
end

if nargin<3,
  flagmoments=0;
end

M_.params(indx) = params;
% [A(oo_.dr.order_var,oo_.dr.order_var),B(oo_.dr.order_var,:)]=dynare_resolve;
[A,B]=dynare_resolve;
if flagmoments==0,
tau = [A(:); vech(B*B')];
else
GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold);
tau = vech(GAM(mf,mf));
end