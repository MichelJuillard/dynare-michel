function time_series = extended_path(initial_conditions,sample_size)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
 % series of size T  is obtained by solving T perfect foresight models. 
 %    
 % INPUTS
 %  o initial_conditions     [double]    m*nlags array, where m is the number of endogenous variables in the model and
 %                                       nlags is the maximum number of lags.
 %  o sample_size            [integer]   scalar, size of the sample to be simulated. 
 %   
 % OUTPUTS
 %  o time_series            [double]    m*sample_size array, the simulations.
 %    
 % ALGORITHM
 %  
 % SPECIAL REQUIREMENTS
 
 % Copyright (C) 2009 Dynare Team
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
 global M_ oo_ options_ 
 
 % Set default initial conditions.
 if isempty(initial_conditions) 
     initial_conditions = repmat(oo_.steady_state,1,M_.maximum_lag); 
 end 
 
 % Copy sample_size to periods.
 options_.periods = sample_size; 
 
 % Initialize the exogenous variables.
 make_ex_; 
 
 % Initialize the endogenous variables.
 make_y_; 
 
 % Initialize the output array.
 time_series = NaN(M_.endo_nbr,sample_size+1); 
 
 % Set the covariance matrix of the structural innovations.
 variances = diag(M_.Sigma_e); 
 positive_var_indx = find(variances>0); 
 covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx); 
 number_of_structural_innovations = length(covariance_matrix); 
 covariance_matrix_upper_cholesky = chol(covariance_matrix); 
 
 tdx = M_.maximum_lag+1; 
 
 for t=1:sample_size 
     oo_.exo_simul(tdx,positive_var_indx) = exp(randn(1,number_of_structural_innovations)*covariance_matrix_upper_cholesky-.5*variances(positive_var_indx)'); 
     perfect_foresight_simulation; 
     time_series(:,t+1) = oo_.endo_simul(:,tdx); 
     oo_.endo_simul(:,1:end-1) = oo_.endo_simul(:,2:end); 
     oo_.endo_simul(:,end) = oo_.steady_state; 
 end 