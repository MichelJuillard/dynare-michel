function f=calib_obj(M_.Sigma_e,A,ghu1,ghx,ghu,targets,var_weights,iy,nar)
% targets and iy order: 1) variances 2) correlations 
% 3) constraints on M_.Sigma_e itself 4) autocorrelations

% Copyright (C) 2005-2008 Dynare Team
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

  global vx fold options_
  
  oo_.gamma_y = cell(nar+1,1);
%  M_.Sigma_e = M_.Sigma_e'*M_.Sigma_e;
  M_.Sigma_e=diag(M_.Sigma_e);
  nx = size(ghx,2);
  b=ghu1*M_.Sigma_e*ghu1';
  vx = [];
  if isempty(vx)
    vx = lyapunov_symm(A,b,options_.qz_criterium,options_.lyapunov_complex_threshold);
  else
    [vx,status] = bicgstab_(@f_var,b(:),vx(:),1e-8,50,A,nx);
    if status
      vx = lyapunov_symm(A,b,options_.qz_criterium,options_.lyapunov_complex_threshold);
    else
      vx=reshape(vx,nx,nx);
    end
  end
  oo_.gamma_y{1} = ghx*vx*ghx'+ ghu*M_.Sigma_e*ghu';
  f = 0;
  if ~isempty(targets{1})
    e = targets{1}-sqrt(oo_.gamma_y{1}(iy{1}));
    f = e'*(var_weights{1}.*e);
  end

  sy = sqrt(diag(oo_.gamma_y{1}));
  sy = sy *sy';
  if ~isempty(targets{2})
    e = targets{2}-oo_.gamma_y{1}(iy{2})./(sy(iy{2})+1e-10);
    f = f+e'*(var_weights{2}.*e);
  end
  
  if ~isempty(targets{3})
    e = targets{3}-sqrt(M_.Sigma_e(iy{3}));
    f = f+e'*(var_weights{3}.*e);
  end
  
  % autocorrelations
  if nar > 0
    vxy = (A*vx*ghx'+ghu1*M_.Sigma_e*ghu');
    
    oo_.gamma_y{2} = ghx*vxy./(sy+1e-10);
    if ~isempty(targets{4})
      e = targets{4}-oo_.gamma_y{2}(iy{4});
      f = f+e'*(var_weights{4}.*e);
    end
    
    for i=2:nar
      vxy = A*vxy;
      oo_.gamma_y{i+1} = ghx*vxy./(sy+1e-10);
      if ~isempty(targets{i+3})
	e = targets{i+3}-oo_.gamma_y{i+1}(iy{i+3});
	f = f+e'*(var_weights{i+3}.*e);
      end
    end
  end
  if isempty(fold) | f < 2*fold
    fold = f;
    vxold = vx;
  end
  % 11/04/02 MJ generalized for correlations, autocorrelations and
  %             constraints on M_.Sigma_e
  % 01/25/03 MJ targets std. dev. instead of variances
  