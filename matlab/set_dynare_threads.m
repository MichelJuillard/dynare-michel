function set_dynare_threads(n)
% This function sets the number of threads used by some MEX files when compiled
% with OpenMP support, i.e with --enable-openmp is given to configure.
% As of 2010-09-27, only A_times_B_kronecker_C and
% sparse_hessian_times_B_kronecker_C support this.
%
% INPUTS 
%  o n    [integer]   scalar specifying the number of threads to be used.    

% Copyright (C) 2009-2010 Dynare Team
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

setenv('DYNARE_NUM_THREADS',int2str(n));
