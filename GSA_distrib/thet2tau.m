function tau = thet2tau(params, indx)
global M_ oo_

if nargin==1,
    indx = [1:M_.param_nbr];
end

M_.params(indx) = params;
% [A(oo_.dr.order_var,oo_.dr.order_var),B(oo_.dr.order_var,:)]=dynare_resolve;
[A,B]=dynare_resolve;
tau = [A(:); vech(B*B')];