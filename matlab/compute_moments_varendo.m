function oo_ = compute_moments_varendo(options_,M_,oo_,var_list_)
% Computes the second order moments (autocorrelation function, covariance
% matrix and variance decomposition) for all the endogenous variables selected in
% var_list_. The results are saved in oo_
%  
% INPUTS:
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   oo_             [structure]    Dynare structure (results).
%   var_list_       [string]       Array of string with endogenous variable names.
%    
% OUTPUTS
%   oo_             [structure]    Dynare structure (results).
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008 Dynare Team
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
    NumberOfEndogenousVariables = rows(var_list_);
    NumberOfExogenousVariables = M_.exo_nbr;
    list_of_exogenous_variables = M_.exo_names;
    NumberOfLags = options_.ar;
    % COVARIANCE MATRIX.
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfEndogenousVariables
            oo_ = posterior_analysis('variance',var_list_(i,:),var_list_(j,:),[],options_,M_,oo_);
        end
    end
    % CORRELATION FUNCTION.
    for h=NumberOfLags:-1:1
        for i=1:NumberOfEndogenousVariables
            for j=1:NumberOfEndogenousVariables
                oo_ = posterior_analysis('correlation',var_list_(i,:),var_list_(j,:),h,options_,M_,oo_);
            end
        end
    end
    % VARIANCE DECOMPOSITION.
    for i=1:NumberOfEndogenousVariables
        for j=i:NumberOfExogenousVariables
            oo_ = posterior_analysis('decomposition',var_list_(i,:),M_.exo_names(j,:),[],options_,M_,oo_);
        end
    end