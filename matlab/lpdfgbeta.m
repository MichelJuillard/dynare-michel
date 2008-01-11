function ldens = lpdfgbeta(x,a,b,aa,bb);

% function ldens = lpdfgbeta(x,a,b,aa,bb);
% log (generalized) BETA PDF 
%
% INPUTS
%    x:      density evatuated at x
%    a:      BETA distribution paramerer 
%    b:      BETA distribution paramerer 
%    aa:     lower bound 
%    bb:     upper bound 

% OUTPUTS
%    ldens:  the log (generalized) BETA PDF
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

ldens = -betaln(a,b) + (a-1)*log(x-aa) + (b-1)*log(bb-x) - (a+b-1)*log(bb-aa);

%gammaln(a+b) - gammaln(a) - gammaln(b)
%betaln(a,b)
%pause
% 02/16/04 SA Interval [aa,bb] is the support of the probability density function. 
