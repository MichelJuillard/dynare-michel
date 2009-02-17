function pdraw = prior_draw(init, cc, npar, prior_structure)
% This function generate one draw from the joint prior distribution.
% 
% INPUTS 
%   o init             [integer]    scalar equal to 1 (first call) or 0.
%   o cc               [double]     two columns matrix (same as in
%                                   metropolis.m), constraints over the
%                                   parameter space (upper and lower bounds).
%   o npar             [integer]    number of estimated parameters.
%   o prior_structure  [structure]  Describes the prior distribution [bayestopt_]
%
% OUTPUTS 
%   o pdraw            [double]     npar*1 vector, draws from the joint prior density.
%
%
% SPECIAL REQUIREMENTS
%   none
%
% NOTE Input arguments 1,2,3 and 4 are only needed for initialization.

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

persistent number_of_estimated_parameters prior_shape prior_mean prior_standard_deviation a b p1 p2 p3 p4
persistent bounds nobounds_flag
persistent uniform_index gaussian_index gamma_index beta_index inverse_gamma_1_index inverse_gamma_2_index

if nargin>0 && init
    number_of_estimated_parameters = npar;
    prior_shape = prior_structure.pshape;
    prior_mean = prior_structure.pmean;
    prior_standard_deviation  = prior_structure.pstdev;
    p1 = prior_structure.p1;
    p2 = prior_structure.p2;
    p3 = prior_structure.p3;
    p4 = prior_structure.p4;
    a = NaN(number_of_estimated_parameters,1);
    b = NaN(number_of_estimated_parameters,1); 
    if ~isempty(cc)
        bounds = cc;
        nobounds_flag = 0;
    else
        nobounds_flag = 1;
    end
    beta_index = find(prior_shape==1);
    gamma_index = find(prior_shape==2);
    gaussian_index = find(prior_shape==3);
    inverse_gamma_1_index = find(prior_shape==4);
    uniform_index = find(prior_shape==5);
    inverse_gamma_2_index = find(prior_shape==6);
    % Set parameters for the beta prior
    mu = (p1(beta_index)-p3(beta_index))./(p4(beta_index)-p3(beta_index));
    stdd = p2(beta_index)./(p4(beta_index)-p3(beta_index));
    a(beta_index) = (1-mu).*mu.^2./stdd.^2 - mu;
    b(beta_index) = a(beta_index).*(1./mu - 1);
    % Set parameters for the gamma prior
    mu = p1(gamma_index)-p3(gamma_index);
    b(gamma_index) = p2(gamma_index).^2./mu;
    a(gamma_index) = mu./b(gamma_index);
    % Initialization of the vector of prior draws.
    pdraw = zeros(number_of_estimated_parameters,1);
    return
end

if nobounds_flag
    % Uniform draws.
    if ~isempty(uniform_index)
        pdraw(uniform_index) = rand(length(uniform_index),1).*(p4(uniform_index)-p3(uniform_index)) + p3(uniform_index);  
    end
    % Gaussian draws.
    if ~isempty(gaussian_index)
        pdraw(gaussian_index) = randn(length(gaussian_index),1).*prior_standard_deviation(gaussian_index) + prior_mean(gaussian_index);
    end
    % Gamma draws.
    if ~isempty(gamma_index)
        pdraw(gamma_index) = gamrnd(a(gamma_index),b(gamma_index))+p3(gamma_index);
    end
    % Beta draws.
    if ~isempty(beta_index)
        pdraw(beta_index) = (p4(beta_index)-p3(beta_index)).*betarnd(a(beta_index),b(beta_index))+p3(beta_index);
    end
    % Inverted gamma (type 1) draws.
    if ~isempty(inverse_gamma_1_index)
        pdraw(inverse_gamma_1_index) = ...
            sqrt(1./gamrnd(p2(inverse_gamma_1_index)/2,2./p1(inverse_gamma_1_index)));
    end
    % Inverted gamma (type 2) draws.
    if ~isempty(inverse_gamma_2_index)
        pdraw(inverse_gamma_2_index) = ...
            1./gamrnd(p2(inverse_gamma_2_index)/2,2./p1(inverse_gamma_2_index));
    end
else
    if ~isempty(uniform_index)
        mm = length(uniform_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(uniform_index(id)) = rand(mm,1).*(p4(uniform_index(id))-p3(uniform_index(id))) + p3(uniform_index(id));
            id = ID( find( (pdraw(uniform_index(id))<bounds(uniform_index(id),1)) || ...
                           (pdraw(uniform_index(id))>bounds(uniform_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
    if ~isempty(gaussian_index)
        mm = length(gaussian_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(gaussian_index(id)) = randn(mm,1).*prior_standard_deviation(gaussian_index(id)) + prior_mean(gaussian_index(id));
            id = ID( find( (pdraw(gaussian_index(id))<bounds(gaussian_index(id),1)) || ...
                           (pdraw(gaussian_index(id))>bounds(gaussian_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
    if ~isempty(gamma_index)
        mm = length(gamma_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(gamma_index(id)) = gamrnd(a(gamma_index(id)),b(gamma_index(id)))+p3(gamma_index(id));
            id = ID( find( (pdraw(gamma_index(id))<bounds(gamma_index(id),1)) || ...
                           (pdraw(gamma_index(id))>bounds(gamma_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
    if ~isempty(beta_index)
        mm = length(beta_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(beta_index(id)) = (p4(beta_index(id))-p3(beta_index(id))).*betarnd(a(beta_index(id)),b(beta_index(id)))+p3(beta_index(id));
            id = ID( find( (pdraw(beta_index(id))<bounds(beta_index(id),1)) || ...
                           (pdraw(beta_index(id))>bounds(beta_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
    if ~isempty(inverse_gamma_1_index)
        mm = length(inverse_gamma_1_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(inverse_gamma_1_index(id)) = sqrt(1./gamrnd(p2(inverse_gamma_1_index(id))/2,2./p1(inverse_gamma_1_index(id))));
            id = ID( find( (pdraw(inverse_gamma_1_index(id))<bounds(inverse_gamma_1_index(id),1)) || ...
                           (pdraw(inverse_gamma_1_index(id))>bounds(inverse_gamma_1_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
    if ~isempty(inverse_gamma_2_index)
        mm = length(inverse_gamma_2_index);
        id = 1:mm;
        ID = 1:mm;
        while mm
            pdraw(inverse_gamma_2_index(id)) = 1./gamrnd(p2(inverse_gamma_2_index(id))/2,2./p1(inverse_gamma_2_index(id)));
            id = ID( find( (pdraw(inverse_gamma_2_index(id))<bounds(inverse_gamma_2_index(id),1)) || ...
                           (pdraw(inverse_gamma_2_index(id))>bounds(inverse_gamma_2_index(id),2)) ) ) ;
            mm = length(id);
        end
    end
end