function pdraw = prior_draw(init,  prior_structure)
% This function generate one draw from the joint prior distribution.
% 
% INPUTS 
%   o init             [integer]    scalar equal to 1 (first call) or 0.
%   o prior_structure  [structure]  Describes the prior distribution [bayestopt_]
%
% OUTPUTS 
%   o pdraw            [double]     1*npar vector, draws from the joint prior density.
%
%
% SPECIAL REQUIREMENTS
%   none
%
% NOTE 1. Input arguments 1 an 2 are only needed for initialization.
% NOTE 2. A given draw from the joint prior distribution does not satisfy BK conditions a priori.

% Copyright (C) 2006-2009 Dynare Team
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

persistent p6 p7 p3 p4
persistent uniform_index gaussian_index gamma_index beta_index inverse_gamma_1_index inverse_gamma_2_index
persistent uniform_draws gaussian_draws gamma_draws beta_draws inverse_gamma_1_draws inverse_gamma_2_draws


if nargin>0 && init
    p6 = prior_structure.p6;
    p7 = prior_structure.p7;
    p3 = prior_structure.p3;
    p4 = prior_structure.p4;
    prior_shape = prior_structure.pshape ;
    number_of_estimated_parameters = length(p6);
    beta_index = find(prior_shape==1);
    beta_draws = 1;
    if isempty(beta_index)
        beta_draws = 0;
    end
    gamma_index = find(prior_shape==2);
    gamma_draws = 1;
    if isempty(gamma_index)
        gamma_draws = 0;
    end
    gaussian_index = find(prior_shape==3);
    gaussian_draws = 1;
    if isempty(gaussian_index)
        gaussian_draws = 0;
    end
    inverse_gamma_1_index = find(prior_shape==4);
    inverse_gamma_1_draws = 1;
    if isempty(inverse_gamma_1_index)
        inverse_gamma_1_draws = 0;
    end
    uniform_index = find(prior_shape==5);
    uniform_draws = 1;
    if isempty(uniform_index)
        uniform_draws = 0;
    end
    inverse_gamma_2_index = find(prior_shape==6);
    inverse_gamma_2_draws = 1;
    if isempty(inverse_gamma_2_index)
        inverse_gamma_2_draws = 0;
    end
    pdraw = zeros(number_of_estimated_parameters,1);
    return
end

if uniform_draws
    pdraw(uniform_index) = rand(length(uniform_index),1).*(p4(uniform_index)-p3(uniform_index)) + p3(uniform_index);  
end

if gaussian_draws
    pdraw(gaussian_index) = randn(length(gaussian_index),1).*p7(gaussian_index) + p6(gaussian_index);
end

if gamma_draws
    pdraw(gamma_index) = gamrnd(p6(gamma_index),p7(gamma_index))+p3(gamma_index);
end

if beta_draws
    pdraw(beta_index) = (p4(beta_index)-p3(beta_index)).*betarnd(p6(beta_index),p7(beta_index))+p3(beta_index);
end

if inverse_gamma_1_draws
    pdraw(inverse_gamma_1_index) = ...
        sqrt(1./gamrnd(p7(inverse_gamma_1_index)/2,2./p6(inverse_gamma_1_index)))+p3(inverse_gamma_1_index);
end

if inverse_gamma_2_draws
    pdraw(inverse_gamma_2_index) = ...
        1./gamrnd(p7(inverse_gamma_2_index)/2,2./p6(inverse_gamma_2_index))+p3(inverse_gamma_2_index);
end