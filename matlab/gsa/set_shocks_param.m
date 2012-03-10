function set_shocks_param(xparam1)

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

global estim_params_ M_
  
  nvx = estim_params_.nvx;
  ncx = estim_params_.ncx;
  np = estim_params_.np;
  Sigma_e = M_.Sigma_e;
  offset = 0;
  if nvx
    offset = offset + nvx;
    var_exo = estim_params_.var_exo;
    for i=1:nvx
      k = var_exo(i,1);
      Sigma_e(k,k) = xparam1(i)^2;
    end
  end
  
  if ncx
    offset = offset + estim_params_.nvn;
    corrx = estim_params_.corrx;
    for i=1:ncx
      k1 = corrx(i,1);
      k2 = corrx(i,2);
      Sigma_e(k1,k2) = xparam1(i+offset)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
      Sigma_e(k2,k1) = Sigma_e_(k1,k2);
    end
  end
  
  
  M_.Sigma_e = Sigma_e;