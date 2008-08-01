function D = sparse_hessian_times_B_kronecker_C(A,B,C)
%function D = sparse_hessian_times_B_kronecker_C(A,B,C)
% Computes A * kron(B,C) where A is a sparse matrix.
%
% INPUTS
%   A  [double] mA*nA matrix.
%   B  [double] mB*nB matrix.
%   C  [double] mC*nC matrix.
%  
% OUTPUTS
%   D  [double] mA*(nC*nB) or mA*(nB*nB) matrix.
%  
% ALGORITHM
%   none.    
%
% SPECIAL REQUIREMENTS
%   none.

% Copyright (C) 1996-2008 Dynare Team
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

switch nargin
  case 3
    D = A_times_B_kronecker_C(A,B,C);
  case 2
    D = A_times_B_kronecker_C(A,B,B);
  otherwise
    error('Two or Three input arguments required!')
end
