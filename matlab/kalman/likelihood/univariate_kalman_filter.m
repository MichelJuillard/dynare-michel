function [LIK, lik,a,P] = univariate_kalman_filter(data_index,number_of_observations,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods)
% Computes the likelihood of a stationnary state space model (univariate approach).
%
% INPUTS
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]    pp*smpl matrix of data.
%    start                        [integer]   scalar, index of the first observation (column of Y).
%    last                         [integer]   scalar, index of the last observation (column of Y).    
%    a                            [double]    mm*1 vector, initial level of the state vector.
%    P                            [double]    mm*mm matrix, covariance matrix of the initial state vector.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).    
%    presample                    [integer]   scalar, presampling if strictly positive.
%    T                            [double]    mm*mm transition matrix of the state equation.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.    
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    H                            [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    Z                            [integer]   pp*1 vector of indices for the observed variables.    
%    mm                           [integer]   scalar, dimension of the state vector.
%    pp                           [integer]   scalar, number of observed variables.
%    rr                           [integer]   scalar, number of structural innovations.
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%    a          [double]    mm*1 vector, estimated level of the states.
%    P          [double]    mm*mm matrix, covariance matrix of the states.    
%
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2011 Dynare Team
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

if nargin<20 || isempty(Zflag)% Set default value for Zflag ==> Z is a vector of indices.
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<21
    diffuse_periods = 0;
end
    
% Get sample size.
smpl = last-start+1;

% Initialize some variables.
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldP = Inf;
l2pi = log(2*pi);
notsteady = 1;

oldK = Inf;
K = NaN(mm,pp);

while notsteady && t<=last
    s = t-start+1;
    d_index = data_index{t};
    if Zflag
        z = Z(d_index,:);
    else
        z = Z(d_index);
    end
    oldP = P(:);
    for i=1:rows(z)
        if Zflag
            prediction_error = Y(d_index(i),t) - z(i,:)*a;
            Fi = z(i,:)*P*z(i,:)' + H(d_index(i),d_index(i));
        else
            prediction_error = Y(d_index(i),t) - a(z(i));
            Fi = P(z(i),z(i)) + H(d_index(i),d_index(i));
        end
        if Fi>kalman_tol
            if Zflag
                Ki =  P*z(i,:)'/Fi;
            else
                Ki = P(:,z(i))/Fi;
            end
            if t>no_more_missing_observations
                K(:,i) = Ki;
            end
            a = a + Ki*prediction_error;
            P = P - (Fi*Ki)*Ki';
            lik(s) = lik(s) + log(Fi) + prediction_error*prediction_error/Fi + l2pi;
        end
    end
    a = T*a;
    P = T*P*transpose(T) + QQ;
    if t>=no_more_missing_observations
        notsteady = max(abs(K(:)-oldK))>riccati_tol;
        oldK = K(:);
    end
    t = t+1;
end

% Divide by two.
lik(1:s) = .5*lik(1:s);

% Call steady state univariate kalman filter if needed.
if t<last
    [tmp, lik(s+1:end)] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag);
end

% Compute minus the log-likelihood.
if presample
    if presample>=diffuse_periods
        lik = lik(1+(presample-diffuse_periods):end);
    end
end
LIK = sum(lik);