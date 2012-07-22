function [x,c] = linsolve(A,B,opts)
% (very imperfect) Clone of Matlab's linsolve.

% Copyright (C) 2010-2011 Dynare Team
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

    c = [];
    x = [];
    if nargin == 3
        if isfield(opts,'TRANSA')
            A = A';
        end
    end
    if nargout == 2
        c = rcond(A);
    end
    
    x = A\B;
    
