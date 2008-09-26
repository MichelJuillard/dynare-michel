function result = check
% function result = check
% checks determinacy conditions by computing the eigenvalues
%
% INPUTS
%    none
%    
% OUTPUTS
%    result
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2008 Dynare Team
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

global M_ options_ oo_
global it_
  
  temp_options = options_;
  tempex = oo_.exo_simul;
  if ~options_.initval_file & M_.exo_nbr > 1
    oo_.exo_simul = ones(M_.maximum_lead+M_.maximum_lag+1,1)*oo_.exo_steady_state';
  end
  
  options_ = set_default_option(options_,'noprint',0);
  options_ = set_default_option(options_,'linear',0);  
  options_ = set_default_option(options_,'steadystate_flag',0);
  options_.order = 1;
  if exist([M_.fname '_steadystate.m'])
    options_.steadystate_flag = 1;
  end

  [dr, info] = resol(oo_.steady_state,1);
  
  oo_.dr = dr;
  
  if info(1) ~= 0 & info(1) ~= 3 & info(1) ~= 4
    print_info(info);
  end  

  oo_.exo_simul = tempex;
  
  eigenvalues_ = dr.eigval;
  if(options_.model_mode == 1)
    nyf = dr.nyf;      
  else
    nyf = nnz(dr.kstate(:,2)>M_.maximum_lag+1);
  end;
  [m_lambda,i]=sort(abs(eigenvalues_));
  n_explod = nnz(abs(eigenvalues_) > options_.qz_criterium);
  
  if options_.noprint == 0
    disp(' ')
    disp('EIGENVALUES:')
    disp(sprintf('%16s %16s %16s\n','Modulus','Real','Imaginary'))
    z=[m_lambda real(eigenvalues_(i)) imag(eigenvalues_(i))]';
    disp(sprintf('%16.4g %16.4g %16.4g\n',z))
    options_ = set_default_option(options_,'qz_criterium',1.000001);
    disp(sprintf('\nThere are %d eigenvalue(s) larger than 1 in modulus ', ...
		 n_explod));
    disp(sprintf('for %d forward-looking variable(s)',nyf));
    disp(' ')
    if dr.rank == nyf & nyf == n_explod
      disp('The rank condition is verified.')
    else
      disp('The rank conditions ISN''T verified!')
    end
    disp(' ')
  end
  
  % keep lambda_ for backward compatibility
  lambda_ = eigenvalues_;

  options_ = temp_options;
  
  % 2/9/99 MJ: line 15, added test for absence of exogenous variable.
  % 8/27/2000 MJ: change JACOB call. Added ...,1 to cumsum()
  % 6/24/01 MJ: added count of abs(eigenvalues) > 1
  % 2/21/02 MJ: count eigenvalues > 1[+1e-5]
  % 01/22/03 MJ: warning(warning_state) needs parentheses for Matlab 6.5
  %              name conflicts with parameters
  % 05/21/03 MJ: replace computation by dr1.m and add rank check
  % 06/05/03 MJ: corrected bug when M_.maximum_lag > 0
  




