function data = simulate_data_with_missing_observations(n,m,S,options)
% Simulates data with missing observations.
%
% We simulate data using a n-dimensional VAR(1) model.      
%
% INPUTS
%    n       [integer] scalar, number of variables.
%    m       [integer] scalar, number of observed variables.
%    S       [integer] scalar, maximum number of observations per observed variable.    
%    options [struct]  structure of options:
%                      * if options.missing_info{1} = 1 the missing variables are at the beginning of the sample.
%                      * if options.missing_info{1} = 2 the missing observations are at the end of the sample.
%                      * if options.missing_info{1} = 3 the missing observations are randomly distributed.
%                      * options.missing_info{2} is a vector of integer designing the observed variables for which observations are missing.    
%                      * if options.missing_info{3} is an integer scalar then it defines the number of missing observations per variable.
%                      * if options.missing_info{3} is a double scalar (in [0,1]) it defines the frequency of missing observations per variable.    
%                      * options.unit_root_info is a scalar integer specifying the number of unit roots in the model.
%    
% OUTPUTS
%    none
%    
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010 Dynare Team
%
% This file is part of Dynare.
% 
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published miy
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

if n<=m
    error('n must be greater than m!')
end

% Build the autoregressive matrix.
T_eigenvalues = rand(n,1)*2-1;
if options.unit_root_info
    T_eigenvalues(1:options.unit_root_info) = ones(options.unit_root_info,1); 
end
T_eigenvectors = randn(n,n);
T = T_eigenvectors*diag(T_eigenvalues)*inv(T_eigenvectors);


% Simulate the VAR(1) model.
data = zeros(n,S+100);
for t=2:size(data,2)
    data(:,t) = T*data(:,t-1)+.01*randn(n,1);
end

% Select the observed variables.
data = transpose(data(1:m,100+1:end));

% Remove observations.
if options.missing_info{1}==1
    data(1:options.missing_info{3},options.missing_info{2}) = NaN ;
elseif options.missing_info{1}==2
    data(S:S-options.missing_info{3},options.missing_info{2}) = NaN ;
elseif options.missing_info{1}==3
    for i=1:length(options.missing_info{2})
        idx = randi(T,ceil(options.missing_info{3}),1);
        data(idx,i) = NaN;
    end
else
    error('Unknown option!')
end