% Copyright (C) 2001 Michel Juillard
%
function y=ff1_(x)
global it_ M_ oo_

n1 = size(x,1) - M_.exo_nbr;
oo_.exo_simul(it_+M_.maximum_lag-M_.maximum_lag,:) = x(n1+1:end)';
fh = str2func([M_.fname '_static']);
y=feval(fh,x(1:n1),oo_.exo_simul);



