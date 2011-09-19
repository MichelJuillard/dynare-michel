function [LIK, lik, a, Pstar, llik] = univariate_kalman_filter_d(data_index, number_of_observations, no_more_missing_observations, ...
                                                              Y, start, last, a, Pinf, Pstar, kalman_tol, riccati_tol, presample, ...
                                                              T,R,Q,H,Z,mm,pp,rr)
% Computes the diffuse likelihood of a stationnary state space model (univariate approach).
%
% INPUTS
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%    Y                            [double]      pp*smpl matrix of (detrended) data, where pp is the number of observed variables.
%    start                        [integer]     scalar, first observation.
%    last                         [integer]     scalar, last observation.
%    a                            [double]      mm*1 vector, levels of the state variables.
%    Pinf                         [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    Pstar                        [double]      mm*mm matrix used to initialize the covariance matrix of the state vector.
%    kalman_tol                   [double]      scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]      scalar, tolerance parameter (riccati iteration).
%    presample                    [integer]     scalar, presampling if strictly positive.    
%    T                            [double]      mm*mm matrix, transition matrix in  the state equations.
%    R                            [double]      mm*rr matrix relating the structural innovations to the state vector.
%    Q                            [double]      rr*rr covariance matrix of the structural innovations.
%    H                            [double]      pp*pp covariance matrix of the measurement errors (if H is equal to zero (scalar) there is no measurement error). 
%    Z                            [double]      pp*mm matrix, selection matrix or pp linear independant combinations of the state vector.
%    mm                           [integer]     scalar, number of state variables.
%    pp                           [integer]     scalar, number of observed variables.
%    rr                           [integer]     scalar, number of structural innovations.    
%
% OUTPUTS
%    dLIK         [double]      scalar, minus loglikelihood
%    dlik         [double]      d*1 vector, log density of each vector of observations (where d is the number of periods of the diffuse filter).
%    a            [double]      mm*1 vector, current estimate of the state vector.
%    Pstar        [double]      mm*mm matrix, covariance matrix of the state vector.    
%    llik         [double]      d*pp matrix used by DsgeLikelihood_hh.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "llik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2011 Dynare Team
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
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
dlik = zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
oldK = Inf;
llik = NaN(smpl,pp);

newRank = rank(Pinf,kalman_tol);
l2pi = log(2*pi);

while newRank && (t<=last)
    s = t-start+1;
    d_index = data_index{t};
    for i=1:length(d_index)
        Zi = Z(d_index(i),:);
        prediction_error = Y(d_index(i),t) - Zi*a;
        Fstar = Zi*Pstar*Zi' + H(d_index(i),d_index(i));
        Finf  = Zi*Pinf*Zi';
        Kstar = Pstar*Zi';
        if Finf>kalman_tol && newRank
            Kinf   = Pinf*Zi';
            Kinf_Finf = Kinf/Finf;
            a         = a + Kinf_Finf*prediction_error;
            Pstar     = Pstar + Kinf*(Kinf_Finf'*(Fstar/Finf)) - Kstar*Kinf_Finf' - Kinf_Finf*Kstar';
            Pinf      = Pinf - Kinf*Kinf_Finf';
            llik(s,d_index(i)) = log(Finf) + l2pi;
            dlik(s) = dlik(s) + llik(s,d_index(i));
        elseif Fstar>kalman_tol
            llik(s,d_index(i)) = log(Fstar) + prediction_error*prediction_error/Fstar + l2pi;
            dlik(s) = dlik(s) + llik(s,d_index(i));
            a = a+Kstar*(prediction_error/Fstar);
            Pstar = Pstar-Kstar*(Kstar'/Fstar);
        end
    end
    if newRank
        oldRank = rank(Pinf,kalman_tol);
    else
        oldRank = 0;
    end
    a     = T*a;
    Pstar = T*Pstar*T'+QQ;
    Pinf  = T*Pinf*T';
    if newRank
        newRank = rank(Pinf,kalman_tol);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')   
    end
    t = t+1;
end

if (t>last)
    error(['univariate_diffuse_kalman_filter:: There isn''t enough information to estimate the initial conditions of the nonstationary variables']);
    LIK = NaN;
    return
end

% Divide by two.
dlik = .5*dlik(1:s);
llik = .5*llik(1:s,:);

if presample 
    if presample>=length(dlik)
        dLIK = 0;
        dlik = [];
        llik = [];
    else
        dlik = dlik(1+presample:end);
        llik = llik(1+presample:end,:);
        dLIK = sum(dlik);% Minus the log-likelihood (for the initial periods).
    end
else
    dLIK = sum(dlik);
end