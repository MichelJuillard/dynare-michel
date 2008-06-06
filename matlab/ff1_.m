function y=ff1_(x)

% function y=ff1_(x)
% splits the input argument x into endogenous and exogenous variables and calls the 'static' function
%
% INPUTS
%    x:          argument splitted between endogenous and exogenous
%        
% OUTPUTS
%    y:         'static' function residuals
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2001-2008)
% Gnu Public License.

global it_ M_ oo_

n1 = size(x,1) - M_.exo_nbr;
oo_.exo_simul(it_+M_.maximum_lag-M_.maximum_lag,:) = x(n1+1:end)';
fh = str2func([M_.fname '_static']);
y=feval(fh,x(1:n1),oo_.exo_simul, M_.params);



