function [LIK, lik,a] = univariate_kalman_filter_ss(Y,start,last,a,P,kalman_tol,T,H,Z,pp,Zflag)
% Computes the likelihood of a stationnary state space model (steady state univariate kalman filter).
%
% INPUTS 
%   Y                       [double]    pp*smpl matrix of data.
%   start                   [integer]   scalar, index of the first observation (column of Y).
%   last                    [integer]   scalar, index of the last observation (column of Y).    
%   a                       [double]    mm*1 vector, initial level of the state vector.
%   P                       [double]    mm*mm matrix, covariance matrix of the initial state vector.
%   T                       [double]    mm*mm transition matrix of the state equation.
%   Z                       [integer]   pp*1 vector of indices for the observed variables.    
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
l2pi = log(2*pi);

% Steady state kalman filter.
while t<=last
    s  = t-start+1;
    PP = P;
    for i=1:pp
        if Zflag
            prediction_error = Y(i,t) - Z(i,:)*a;
            Fi = Z(i,:)*PP*Z(i,:)' + H(i,i);
        else
            prediction_error = Y(i,t) - a(Z(i));
            Fi = PP(Z(i),Z(i)) + H(i,i);
        end
        if Fi>kalman_tol
            if Zflag
                Ki = PP*Z(i,:)'/Fi;
            else
                Ki = PP(:,Z(i))/Fi;
            end
            a  = a + Ki*prediction_error;
            PP = PP - (Fi*Ki)*transpose(Ki);
            lik(s) = lik(s) + log(Fi) + prediction_error*prediction_error/Fi + l2pi;
        end
    end
    a = T*a;
    t = t+1;
end

lik = .5*lik;

LIK = sum(lik);