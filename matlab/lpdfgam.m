function  ldens = lpdfgam(x,a,b);

% function ldens = lpdfbeta(x,a,b)
% log GAMMA PDF
%
% INPUTS
%    x:      density evatuated at x
%    a:      GAMMA distribution paramerer 
%    b:      GAMMA distribution paramerer 
%    
% OUTPUTS
%    ldens:  the log GAMMA PDF
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  ldens = -gammaln(a) -a*log(b)+ (a-1)*log(x) -x/b ;

% 10/11/03  MJ adapted from an earlier GAUSS version by F. Schorfeide,
%              translated to MATLAB by R. Wouters  
%              use MATLAB gammaln rather than lngam
