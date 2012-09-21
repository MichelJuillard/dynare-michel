function [g, badg] = numgrad2_(fcn,f0,x,epsilon,scale,varargin)
% function [g badg] = numgrad2(fcn,xvarargin)

% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/numgrad.m

% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2012 Dynare Team
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

fh = NaN;

delta = epsilon*max(abs(x));
n = length(x);
g = zeros(n,1);

badg=0;
goog=1;
g0 = 0;

for i=1:n
    xiold = x(i);
    h = step_length_correction(xiold,scale,i)*delta;
    x(i) = xiold + h;
    [fh,junk1,junk2,cost_flag] = feval(fcn, x, varargin{:});
    if cost_flag
        g0 = (fh - f0)/h;
    else
        x(i) = xiold - h;
        [fh,junk1,junk2,cost_flag] = feval(fcn, x, varargin{:});
        if cost_flag
            g0 = (f0-fh)/h;
        else
            goog = 0;
        end
    end
    if goog && abs(g0)< 1e15
        g(i) = g0;
    else
        disp('bad gradient ------------------------')
        % fprintf('Gradient w.r.t. %3d: %10g\n',i,g0)
        g(i) = 0;
        badg = 1;
    end
    x(i) = xiold;
end