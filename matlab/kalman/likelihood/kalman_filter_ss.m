function [LIK, lik, a] = kalman_filter_ss(Y,start,last,a,T,K,iF,dF,Z,pp,Zflag)
% Computes the likelihood of a stationnary state space model (steady state kalman filter).
%
% INPUTS 
%   Y                       [double]    pp*smpl matrix of data.
%   start                   [integer]   scalar, index of the first observation (column of Y).
%   last                    [integer]   scalar, index of the last observation (column of Y).    
%   a                       [double]    mm*1 vector, initial level of the state vector.
%   P                       [double]    mm*mm matrix, covariance matrix of the initial state vector.
%   T                       [double]    mm*mm transition matrix of the state equation.
%   K                       [double]    mm*pp matrix, steady state kalman gain.    
%   iF                      [double]    pp*pp matrix, inverse of the steady state covariance matrix of the predicted errors.
%   dF                      [double]    scalar, determinant of the steady state covariance matrix of the predicted errors. 
%   Z                       [integer]   pp*1 vector of indices for the observed variables, if Zflag=0.    
%   pp                      [integer]   scalar, number of observed variables.
%
%
% OUTPUTS 
%    LIK        [double]    scalar, minus log likelihood.    
%    lik        [double]    (last-start+1)*1 vector, density of each observation.
%    a          [double]    mm*1 vector, estimate of the state vector.     
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2011 Dynare Team
% stephane DOT adjemian AT ens DOT fr
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

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.

while t <= last
    if Zflag
        v = Y(:,t)-Z*a;
    else
        v = Y(:,t)-a(Z);
    end
    a = T*(a+K*v);
    lik(t-start+1) = transpose(v)*iF*v;
    t = t+1;
end

% Adding constant determinant of F (prediction error covariance matrix)
lik = lik + log(dF);

% Add log-likelihhod constants and divide by two
lik = .5*(lik + pp*log(2*pi));

% Sum the observation's densities (minus the likelihood) 
LIK = sum(lik);