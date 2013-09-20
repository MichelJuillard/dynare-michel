function [g, badg] = numgrad2(fcn,f0,x,epsilon,varargin)
% function [g badg] = numgrad2(fcn,xvarargin)

% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/numgrad.m

% Copyright (C) 1993-2007 Christopher Sims
% Copyright (C) 2008-2013 Dynare Team
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

delta = epsilon;
n=length(x);
tvec=delta*eye(n);
g=zeros(n,1);

badg=0;
scale=1;
g0 = 0;

right_derivative = 1;

for i=1:n
    if size(x,1)>size(x,2)
        tvecv=tvec(i,:);
    else
        tvecv=tvec(:,i);
    end
    if right_derivative
        [fh,junk1,junk2,cost_flag] = feval(fcn, x+scale*transpose(tvecv), varargin{:});
        g0 = (fh - f0) / (scale*delta);
    else
        [fh,junk1,junk2,cost_flag] = feval(fcn, x-scale*transpose(tvecv), varargin{:});
        g0 = (f0-fh) / (scale*delta);
    end
    if abs(g0)< 1e15
        g(i) = g0;
    else
        g(i) = 0;
        badg = 1;
    end
end