function [param,sigma] = simulated_moments_estimation(xparam,dataset,options)
% Performs estimation by Simulated Moments Method.
%
% INPUTS:
%  xparam          [double]  p*1 vector of initial values for the estimated parameters. 
%  dataset         [      ]  Structure describing the data set.
%  options         [      ]  Structure defining options for SMM.
%
% OUTPUTS: 
%  param           [double]  p*1 vector of point estimates for the parameters.
%  sigma           [double]  p*p covariance matrix of the SMM estimates.
%
% SPECIAL REQUIREMENTS
%  The user has to provide a file where the moment conditions are defined.

% Copyright (C) 2010 Dynare Team
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
% Copyright (C) 2010 Dynare Team
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

global M_
    
% Load the dataset.
eval(dataset.name);
dataset.data = [];
for v = 1:dataset.number_of_observed_variables
    eval(['dataset.data = [ dataset.data , ' dataset.variables(v,:) ' ];'])
end
data = dataset.data(dataset.first_observation:dataset.first_observation+dataset.number_of_observations,:);

% Compute sample moments and the weighting matrix.
eval(['[sample_moments,weighting_matrix] = ' M_.fname '_moments;'])
weighting_matrix = inv(weighting_matrix);

% Initialize output.
sigma = [];
param = [];

% Set options for csminwel. 
H0 = 1e-4*eye(options.estimated_parameters.nb);
ct = 1e-4;
it = 1000;
vb = 2;
            
% Call csminwel.
[fval,param,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
                csminwel('smm_objective',xparam,H0,[],ct,it,2,sample_moments,weighting_matrix,options);