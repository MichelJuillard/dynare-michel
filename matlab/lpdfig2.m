function ldens = lpdfig2(x,s,nu)

% function ldens = lpdfig2(x,s,nu)
% log INVERSE GAMMA (type 2) 
% X ~ IG2(s,nu)
% X = inv(Z) where Z ~ G(nu/2,2/s) (Gamma distribution) 
%
% INPUTS
%    x:      density evatuated at x
%    s:      shape parameter 
%    nu:     scale parameter 

% OUTPUTS
%    ldens:  the log INVERSE GAMMA density function (type 2)
%        
% SPECIAL REQUIREMENTS
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
% details.
%  
% part of DYNARE, copyright Dynare Team (2004-2008)
% Gnu Public License.


ldens = - gammaln(nu/2) - (nu/2)*log(2/s) - .5*(nu+2)*log(x) -.5*s./x;