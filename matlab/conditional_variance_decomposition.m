function PackedConditionalVarianceDecomposition = conditional_variance_decomposition(StateSpaceModel, Steps, SubsetOfVariables)
% This function computes the conditional variance decomposition of a given state space model
% for a subset of endogenous variables.
% 
% INPUTS 
%   StateSpaceModel     [structure]   Specification of the state space model.
%   Steps               [integer]     1*h vector of dates.
%   SubsetOfVariables   [integer]     1*q vector of indices.
%    
% OUTPUTS 
%   PackedConditionalVarianceDecomposition  [double] n(n+1)/2*p matrix, where p is the number of state innovations and
%                                                    n is equal to length(SubsetOfVariables).    
%
% SPECIAL REQUIREMENTS
%
% [1] The covariance matrix of the state innovations needs to be diagonal.
% [2] In this version, absence of measurement errors is assumed...

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
    ConditionalVariance = zeros(StateSpaceModel.number_of_state_equations,StateSpaceModel.number_of_state_equations);
    ConditionalVariance = repmat(ConditionalVariance,[1 1 length(Steps) StateSpaceModel.number_of_state_innovations]);
    BB = StateSpaceModel.impulse_matrix*transpose(StateSpaceModel.impulse_matrix);
    for h = 1:length(Steps)
        for t = 0:Steps(h)
            for i=1:StateSpaceModel.number_of_state_innovations
                ConditionalVariance(:,:,h,i) = ... 
                    StateSpaceModel.transition_matrix*ConditionalVariance(:,:,h,i)*transpose(StateSpaceModel.transition_matrix) ...
                    +BB*StateSpaceModel.state_innovations_covariance_matrix(i,i);
            end
        end
    end
    ConditionalVariance = ConditionalVariance(SubsetOfVariables,SubsetOfVariables,:,:);
    NumberOfVariables = length(SubsetOfVariables);
    PackedConditionalVarianceDecomposition = zeros(NumberOfVariables*(NumberOfVariables+1)/2,length(Steps),StateSpaceModel.number_of_state_innovations); 
    for i=1:StateSpaceModel.number_of_state_innovations
        for h = 1:length(Steps)
            PackedConditionalVarianceDecomposition(:,h,i) = vech(ConditionalVariance(:,:,h,i));
        end
    end