function [LIK, lik] = univariate_diffuse_kalman_filter(T,R,Q,H,Pinf,Pstar,Y,start,Z,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations)
% Computes the likelihood of a stationnary state space model (univariate approach).
%
% INPUTS
%    T                            [double]    mm*mm transition matrix of the state equation.
%    R                            [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                            [double]    rr*rr covariance matrix of the structural innovations.
%    H                            [double]    pp*1 (zeros(pp,1) if no measurement errors) variances of the measurement errors. 
%    P                            [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                            [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                        [integer]   scalar, likelihood evaluation starts at 'start'.
%    Z                            [double]    pp*mm, selection matrix or pp independant linear combinations.
%    kalman_tol                   [double]    scalar, tolerance parameter (rcond).
%    riccati_tol                  [double]    scalar, tolerance parameter (riccati iteration).
%    data_index                   [cell]      1*smpl cell of column vectors of indices.
%    number_of_observations       [integer]   scalar.
%    no_more_missing_observations [integer]   scalar.
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2008 Dynare Team
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
global options_

[pp ,smpl] = size(Y); 
mm = size(T,1);
a  = zeros(mm,1);
QQ = R*Q*transpose(R);
t  = 0;
lik = zeros(smpl+1,1);
lik(smpl+1) = number_of_observations*log(2*pi);
notsteady = 1;
crit = 1.e-6;
newRank	= rank(Pinf,crit);
icc=0;

while newRank && (t<smpl)
    t = t+1;
    Za = Z(data_index{t},:);
    for i=1:length(data_index{t})
        Zi = Z(data_index{t}(i),:);
        prediction_error = Y(data_index{t}(i),t) - Zi*a;
        Fstar = Zi*Pstar*Zi' + H(i);
        Finf  = Zi*Pinf*Zi';
        Kstar = Pstar*Zi';
        if Finf>kalman_tol && newRank
            icc=icc+1;
            Kinf   = Pinf*Zi';
            a	   = a + Kinf*(prediction_error/Finf);
            Pstar  = Pstar + Kinf*(Kinf'*(Fstar/(Finf*Finf))) - (Kstar*Kinf'+Kinf*Kstar')/Finf;
            Pinf   = Pinf - Kinf*(Kinf'/Finf);
            lik(t) = lik(t) + log(Finf);
            if ~isempty(options_.diffuse_d)  
                newRank = (icc<options_.diffuse_d);  
                if newRank && (any(diag(Za*Pinf*Za')>kalman_tol)==0 & rank(Pinf,crit)==0); 
                    options_.diffuse_d = icc;
                    newRank=0;
                    disp('WARNING: Change in OPTIONS_.DIFFUSE_D in univariate DKF')
                    disp(['new OPTIONS_.DIFFUSE_D = ',int2str(icc)])
                    disp('You may have to reset the optimisation')
                end
            else
                newRank = (any(diag(Za*Pinf*Za')>kalman_tol) | rank(Pinf,crit));                 
                if newRank==0
                    P0=	T*Pinf*T';
                    newRank = (any(diag(Za*P0*Za')>kalman_tol) | rank(P0,crit));
                    if newRank==0
                        options_.diffuse_d = icc;
                    end
                end
            end
        elseif Fstar>kalman_tol
            lik(t) = lik(t) + log(Fstar) + prediction_error*prediction_error/Fstar;
            a = a + Kstar*prediction_error/Fstar;
            Pstar = Pstar - Kstar*Kstar'/Fstar;
        end
    end
    if newRank
        oldRank = rank(Pinf,crit);
    else
        oldRank = 0;
    end
    a 	  = T*a;
    Pstar = T*Pstar*T'+QQ;
    Pinf  = T*Pinf*T';
    if newRank
        newRank = rank(Pinf,crit);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')	
    end
end

if (t==smpl)
  error(['univariate_diffuse_kalman_filter:: There isn''t enough information to estimate the initial conditions of the nonstationary variables']);
end

while notsteady && (t<smpl)
    t = t+1;
    oldP = Pstar;
    for i=1:length(data_index{t})
        Zi = Z(data_index{t}(i),:);
        prediction_error = Y(data_index{t}(i),t) - Zi*a;
        Fi   = Zi*Pstar*Zi' + H(i);
        if Fi > kalman_tol
            Ki	= Pstar*Zi';
            a	   = a + Ki*prediction_error/Fi;
            Pstar  = Pstar - Ki*Ki'/Fi;
            lik(t) = lik(t) + log(Fi) + prediction_error*prediction_error/Fi;
        end
    end	
    a 	  = T*a;
    Pstar = T*Pstar*T' + QQ;
    if t>no_more_missing_observations
        notsteady = max(max(abs(Pstar-oldP)))>riccati_tol;
    end
end

while t < smpl
    t = t+1;
    Pstar = oldP;
    for i=1:pp
        Zi = Z(i,:);
        prediction_error = Y(i,t) - Zi*a;
        Fi   = Zi*Pstar*Zi'+H(i);
        if Fi > crit
            Ki 		= Pstar*Zi';
            a 		= a + Ki*prediction_error/Fi;
            Pstar 	= Pstar - Ki*Ki'/Fi;
            lik(t)    	= lik(t) + log(Fi) + prediction_error*prediction_error/Fi;
        end
    end	
    a = T*a;
end

LIK = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);