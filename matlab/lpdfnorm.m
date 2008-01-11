function  f = lpdfnorm(x,m,s)

% function f = lpdfnorm(x,m,s)
% The log of the normal density function 
%
% INPUTS
%    x:      density evatuated at x
%    m:      mean 
%    s:      standard deviation 

% OUTPUTS
%    f:      the log of the normal density function
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.


if nargin<3, s=1; end
if nargin<2, m=0; end
f = -log(s)-log(2*pi)/2-((x-m)./s).^2/2;

