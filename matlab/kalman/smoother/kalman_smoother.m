function [alphahat,epsilonhat,etahat,atilde,P,aK,PK,decomp] = kalman_smoother(T,R,Q,H,P0,Y,start,mf,kalman_tol,riccati_tol)

% function [alphahat,epsilonhat,etahat,a,aK,PK,decomp] = kalman_smoother(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol)
% Computes the kalman smoother of a stationary state space model.
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    P0                     [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    kalman_tol             [double]    scalar, tolerance parameter (rcond).
%    riccati_tol            [double]    scalar, tolerance parameter (riccati iteration).
%
%             
% OUTPUTS
%    alphahat: smoothed state variables (a_{t|T})
%    etahat:   smoothed shocks
%    atilde:   matrix of updated variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t})

% SPECIAL REQUIREMENTS
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98). 

% Copyright (C) 2004-2010 Dynare Team
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

option_filter_covariance = options_.filter_covariance;
option_filter_decomposition = options_.filter_decomposition;

nk   = options_.nk;
smpl = size(Y,2);                               % Sample size.
mm   = size(T,2);                               % Number of state variables.
pp   = size(Y,1);                               % Maximum number of
                                                % observed variables.
rr   = size(Q,1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl+1);
atilde          = zeros(mm,smpl);
K               = zeros(mm,pp,smpl);
aK              = zeros(nk,mm,smpl+nk);
iF              = zeros(pp,pp,smpl);
P               = zeros(mm,mm,smpl+1);
QQ              = R*Q*R';
QRt             = Q*R';
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
epsilonhat      = zeros(rr,smpl);
r               = zeros(mm,smpl+1);
oldK            = 0;

if option_filter_covariance
    PK          = zeros(nk,mm,mm,smpl+nk);
else
    PK          = [];
end

if option_filter_decomposition
    decomp = zeros(nk,mm,rr,smpl+nk);
else
    decomp = [];
end

P(:,:,1) = P0;

t = 0;
notsteady = 1;
F_singular = 1;
while notsteady & t<smpl
    t  = t+1;
    v(:,t)  = Y(:,t)-a(mf,t);
    F  = P(mf,mf,t) + H;
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            atilde(:,t) = a(:,t);
            a(:,t+1) = T*a(:,t);
            P(:,:,t+1) = T*P(:,:,t)*T'+QQ;
        end
    else
        F_singular = 0;
        iF(:,:,t)     = inv(F);
        K1            = P(:,mf,t)*iF(:,:,t);
        atilde(:,t)   = a(:,t) + K1*v(:,t);
        K(:,:,t)      = T*K1;
        a(:,t+1)      = T*atilde(:,t);
        P(:,:,t+1)    = (T*P(:,:,t)-K(:,:,t)*P(mf,:,t))*T'+QQ;
    end
    aK(1,:,t+1) = a(:,t+1);
    if option_filter_covariance
        Pf = P(:,:,t);
        Pf = T*Pf*T' + QQ;
        PK(1,:,:,t+1) = Pf;
    end        
    for jnk=2:nk,
        aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        if option_filter_covariance
            Pf = T*Pf*T' + QQ;
            PK(jnk,:,:,t+jnk) = Pf;
        end
    end
    notsteady = max(max(abs(K(:,:,t)-oldK))) > riccati_tol;
    oldK = K(:,:,t);
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

if t < smpl
    t0 = t;
    while t < smpl
        t = t+1;
        v(:,t) = Y(:,t)-a(mf,t);
        atilde(:,t)   = a(:,t) + K1*v(:,t);
        a(:,t+1)      = T*atilde(:,t);
        aK(1,:,t+1) = a(:,t+1);
        if option_filter_covariance
            Pf = P(:,:,t);
            Pf = T*Pf*T' + QQ;
            PK(1,:,:,t+1) = Pf;
        end        
        for jnk=2:nk,
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
            if option_filter_covariance
                Pf = T*Pf*T' + QQ;
                PK(jnk,:,:,t+jnk) = Pf;
            end
        end
    end
    
    K= cat(3,K(:,:,1:t0),repmat(K(:,:,t0),[1 1 smpl-t0+1]));
    P  = cat(3,P(:,:,1:t0),repmat(P(:,:,t0),[1 1 smpl-t0+1]));
    iF = cat(3,iF(:,:,1:t0),repmat(iF(:,:,t0),[1 1 smpl-t0+1]));
end    

t = smpl+1;
while t>1
    t = t-1;
    r(:,t)        = T'*r(:,t+1);    
    r(mf,t)       = r(mf,t)+iF(:,:,t)*v(:,t) - K(:,:,t)'*r(:,t+1);
    alphahat(:,t) = a(:,t) + P(:,:,t)*r(:,t);
    etahat(:,t)   = QRt*r(:,t);
end
epsilonhat = Y-alphahat(mf,:);

if option_filter_decomposition
    ZRQinv = inv(QQ(mf,mf));
    for t = 1:smpl
        % calculate eta_tm1t
        eta = QRt(:,mf)*iF(:,:,t)*v(:,t);
        AAA = P(:,mf,t)*ZRQinv*bsxfun(@times,R(mf,:),eta');
        % calculate decomposition
        Ttok = eye(mm,mm); 
        decomp(1,:,:,t+1) = AAA;
        for h = 2:nk
            AAA = T*AAA;
            decomp(h,:,:,t+h) = AAA;
        end
    end
end

if ~option_filter_covariance
    P = [];
end

