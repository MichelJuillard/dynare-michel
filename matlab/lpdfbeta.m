function ldens = lpdfbeta(x,a,b);

% function ldens = lpdfbeta(x,a,b)
% log Beta PDF
%
% INPUTS
%    x:      density evatuated at x
%    a:      Beta distribution paramerer 
%    b:      Beta distribution paramerer 
%    
% OUTPUTS
%    ldens:  the log Beta PDF
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  ldens = gammaln(a+b) - gammaln(a) - gammaln(b) + (a-1)*log(x) + (b-1)*log(1-x);

% 10/11/03 MJ adapted from a GAUSS version by F. Schorfheide, translated
%             to Matlab by R. Wouters.  
%             use Matlab gammaln instead of lngam