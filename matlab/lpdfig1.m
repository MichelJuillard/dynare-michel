function ldens = lpdfig1(x,s,nu)
% function ldens = lpdfig1(x,s,nu)
% log INVERSE GAMMA (type 1) 
% X ~ IG1(s,nu)
% X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
%
% INPUTS
%    x:      density evatuated at x
%    s:      shape parameter 
%    nu:     scale parameter 
%
% OUTPUTS
%    ldens:  the log INVERSE GAMMA density function (type 1)
%        
% SPECIAL REQUIREMENTS
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
% details.

% Copyright (C) 2004-2008 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

ldens = log(2) - gammaln(nu/2) - (nu/2).*log(2/s) - (nu+1)*log(x) - .5*s./(x.^2);
