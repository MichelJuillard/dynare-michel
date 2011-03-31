function irf_shocks_indx=getIrfShocksIndx()

% Copyright (C) 2011 Dynare Team
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

global M_ options_

if (isfield(options_,'irf_shocks')==0)
    irf_shocks_indx = M_.exo_names_orig_ord;
else
    irf_shocks_indx = zeros(1,size(options_.irf_shocks,1));
    for i=1:size(options_.irf_shocks,1)
        irf_shocks_indx(i) = find(M_.exo_names==options_.irf_shocks(i));
    end
end