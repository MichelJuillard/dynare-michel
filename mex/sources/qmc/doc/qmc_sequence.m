% Sobol sequence of quasi random numbers.
%
% This function produces a quasi monte carlo sequences.
%
% INPUT (for initialization)
%  o dimension  [uint32 scalar]   dimension {2,3,...,1111}
%  o flag       [uint32 scalar]   flag {0,1,2,3}
%  o seed       [uint32 scalar]   seed (needed if flag>0) 
%  o transform  [uint32 scalar]   transform {0,1}
%
% OUTPUT (initialization)
%
%  o qmc       [structure]  
%    *** qmc.dimension                     [uint32 scalar]
%    *** qmc.flag                          [uint32 scalar]
%    *** qmc.transform                     [uint32 scalar]
%    *** qmc.iteration                     [uint32 scalar]
%    *** qmc.table_of_direction_numbers    [uint32 matrix]
%    *** qmc.table_common_denominator      [uint32 scalar]
%    *** qmc.last                          [double matrix]
%    *** qmc.seed                          [uint32 scalar]
%
% INPUT ()
%  o qmc                   [structure]
%  o number_of_simulations [uint32 scalar] number_of_simulations*dimension 
%
% OUTPUT ()
%  o sArray [double matrix] 
%  o qmc    [structure]
%
% ALGORITHMS 
%       Algorithm 659, Collected Algorithms From ACM. Published in
%       Transactions on Mathematical Software, Vol. 14, n°1, P.88.

% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare (can be used outside Dynare).
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