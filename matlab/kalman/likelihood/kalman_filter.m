function [LIK, lik, a, P] = kalman_filter(Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods)
% Computes the likelihood of a stationnary state space model.
%
% INPUTS 
%    Y                      [double]    pp*smpl matrix of data.
%    start                  [integer]   scalar, index of the first observation (column of Y).
%    last                   [integer]   scalar, index of the last observation (column of Y).    
%    a                      [double]    mm*1 vector, initial level of the state vector.
%    P                      [double]    mm*mm matrix, covariance matrix of the initial state vector.
%    kalman_tol             [double]    scalar, tolerance parameter (rcond).
%    riccati_tol            [double]    scalar, tolerance parameter (riccati iteration).    
%    presample              [integer]   scalar, presampling if strictly positive.
%    T                      [double]    mm*mm transition matrix of the state equation.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.    
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    Z                      [integer]   pp*1 vector of indices for the observed variables, if Zflag=0, pp*mm matrix if Zflag>0.    
%    mm                     [integer]   scalar, dimension of the state vector.
%    pp                     [integer]   scalar, number of observed variables.
%    rr                     [integer]   scalar, number of structural innovations.
%
% OUTPUTS 
%    LIK        [double]    
%    lik        [double]    (last-start+1)*1 vector, density of observations in each periods.
%    a          [double]    mm*1 vector, estimated level of the states.
%    P          [double]    mm*mm matrix, covariance matrix of the states.    
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

%     
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

% Set defaults.
if nargin<17
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<18
    diffuse_periods = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end
    
% Get sample size.
smpl = last-start+1;

% Initialize some variables.
dF   = 1;
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
lik  = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = 1;

while notsteady && t<=last
    s = t-start+1;
    if Zflag
        v  = Y(:,t)-Z*a;
        F  = Z*P*Z' + H;        
    else
        v  = Y(:,t)-a(Z);
        F  = P(Z,Z) + H;
    end
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            a = T*a;
            P = T*P*transpose(T)+QQ;
        end
    else
        F_singular = 0;
        dF     = det(F);
        iF     = inv(F);
        lik(s) = log(dF)+transpose(v)*iF*v;
        if Zflag
            K = P*Z'*iF;
            P = T*(P-K*Z*P)*transpose(T)+QQ;
        else
            K = P(:,Z)*iF;
            P = T*(P-K*P(Z,:))*transpose(T)+QQ;
        end
        a = T*(a+K*v);
        notsteady = max(abs(K(:)-oldK))>riccati_tol;
        oldK = K(:);
    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Add observation's densities constants and divide by two.
lik(1:s) = .5*(lik(1:s) + pp*log(2*pi));

% Call steady state Kalman filter if needed.
if t<last
    [tmp, lik(s+1:end)] = kalman_filter_ss(Y,t,last,a,T,K,iF,dF,Z,pp,Zflag);
end

% Compute minus the log-likelihood.
if presample
    if presample>=diffuse_periods
        lik = lik(1+(presample-diffuse_periods):end);
    end
end
LIK = sum(lik);