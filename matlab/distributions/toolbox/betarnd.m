function rnd = betarnd(a, b)
% BETARND  Random samples from the Beta distribution
%  RND = betarnd(A,B) returns a random sample from the
%  Beta distribution with parameters A and B (i.e. mean of
%  the distribution is A/(A+B) and variance is
%  A*B/(A+B)^2/(A+B+1) ).

% Copyright (C) 2008 Dynare Team
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

if (nargin ~= 2)
  error('betarnd: you must give two arguments');
end

if (~isscalar(a) || ~isscalar(b))
  error('betarnd: A and B should be scalar parameters');
end

if (a <= 0 || a == Inf || b <= 0 || b == Inf)
  rnd = NaN;
else
  x = gamrnd(a, 1);
  rnd = x/(x+gamrnd(b, 1));
end
