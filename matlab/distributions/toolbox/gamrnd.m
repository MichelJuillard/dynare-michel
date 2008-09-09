function rnd = gamrnd(a,b)
% GAMRND  Random samples from the Gamma distribution
%  RND = gamrnd(A,B) returns a random sample from the
%  Gamma distribution with parameters A and B (i.e. mean of
%  the distribution is A*B and variance is A*B^2).
%
% Algorithm of Bauwens, Lubrano & Richard (page 316)

% Copyright (C) 2006-2008 Dynare Team
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
  error('gamrnd: you must give two arguments');
end

if (~isscalar(a) || ~isscalar(b))
  error('gamrnd: A and B should be scalar parameters');
end

if (a <= 0 || a == Inf || b <= 0 || b == Inf)
  rnd = NaN;
  return
end

if a >30
    z = randn;
    rnd = b*(z+sqrt(4*a-1))^2/4; 
else
    condi = 1;
    while condi
        x = -1;
        while x<0
            u1 = rand;
            y = tan(pi*u1);
            x = y*sqrt(2*a-1)+a-1; 
        end
        u2 = rand;
        if log(u2) <= log(1+y^2)+(a-1)*log(x/(a-1))-y*sqrt(2*a-1);
            break
        end
    end
    rnd = x*b;
end
