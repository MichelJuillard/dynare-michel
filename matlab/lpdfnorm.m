function  f = lpdfnorm(x,m,s)
% function f = lpdfnorm(x,m,s)
% The log of the normal density function 
%
% INPUTS
%    x:      density evatuated at x
%    m:      mean 
%    s:      standard deviation 
%
% OUTPUTS
%    f:      the log of the normal density function
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2008 Dynare Team
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

if nargin<3, s=1; end
if nargin<2, m=0; end
f = -log(s)-log(2*pi)/2-((x-m)./s).^2/2;

