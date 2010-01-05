function [alphahat,etahat,atilde, aK] = DiffuseKalmanSmoother1(T,R,Q,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)

% function [alphahat,etahat,a, aK] = DiffuseKalmanSmoother1(T,R,Q,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% Computes the diffuse kalman smoother without measurement error, in the case of a non-singular var-cov matrix 
%
% INPUTS
%    T:        mm*mm matrix
%    R:        mm*rr matrix
%    Q:        rr*rr matrix
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    trend
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%    mf:       observed variables index in the state vector
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

% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric 


global options_

nk = options_.nk;
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl+1);
atilde          = zeros(mm,smpl);
aK              = zeros(nk,mm,smpl+1);
iF              = zeros(pp,pp,smpl);
Fstar           = zeros(pp,pp,smpl);
iFinf           = zeros(pp,pp,smpl);
K               = zeros(mm,pp,smpl);
L               = zeros(mm,mm,smpl);
Linf            = zeros(mm,mm,smpl);
Kstar           = zeros(mm,pp,smpl);
P               = zeros(mm,mm,smpl+1);
Pstar           = zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
crit            = options_.kalman_tol;
crit1       = 1.e-8;
steady          = smpl;
rr              = size(Q,1);
QQ              = R*Q*transpose(R);
QRt                     = Q*transpose(R);
alphahat        = zeros(mm,smpl);
etahat          = zeros(rr,smpl);
r               = zeros(mm,smpl+1);

Z = zeros(pp,mm);
for i=1:pp;
    Z(i,mf(i)) = 1;
end

t = 0;
while rank(Pinf(:,:,t+1),crit1) & t<smpl
    t = t+1;
    v(:,t)                      = Y(:,t) - a(mf,t) - trend(:,t);
    if rcond(Pinf(mf,mf,t)) < crit
        return          
    end
    iFinf(:,:,t)        = inv(Pinf(mf,mf,t));
    PZI                 = Pinf(:,mf,t)*iFinf(:,:,t);
    atilde(:,t)         = a(:,t) + PZI*v(:,t);
    Kinf(:,:,t)         = T*PZI;
    a(:,t+1)            = T*atilde(:,t);
    aK(1,:,t+1)         = a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) = T^(jnk-1)*a(:,t+1);
    end
    Linf(:,:,t)         = T - Kinf(:,:,t)*Z;
    Fstar(:,:,t)        = Pstar(mf,mf,t);
    Kstar(:,:,t)        = (T*Pstar(:,mf,t)-Kinf(:,:,t)*Fstar(:,:,t))*iFinf(:,:,t);
    Pstar(:,:,t+1)      = T*Pstar(:,:,t)*transpose(T)-T*Pstar(:,mf,t)*transpose(Kinf(:,:,t))-Kinf(:,:,t)*Pinf(mf,mf,t)*transpose(Kstar(:,:,t)) + QQ;
    Pinf(:,:,t+1)       = T*Pinf(:,:,t)*transpose(T)-T*Pinf(:,mf,t)*transpose(Kinf(:,:,t));
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
iFinf = iFinf(:,:,1:d);
Linf  = Linf(:,:,1:d);
Fstar = Fstar(:,:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
notsteady = 1;
while notsteady & t<smpl
    t = t+1;
    v(:,t)      = Y(:,t) - a(mf,t) - trend(:,t);
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    if rcond(P(mf,mf,t)) < crit
        return          
    end    
    iF(:,:,t)   = inv(P(mf,mf,t));
    PZI         = P(:,mf,t)*iF(:,:,t);
    atilde(:,t) = a(:,t) + PZI*v(:,t);
    K(:,:,t)    = T*PZI;
    L(:,:,t)    = T-K(:,:,t)*Z;
    a(:,t+1)    = T*atilde(:,t);
    aK(1,:,t+1)                 = a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk)                 = T^(jnk-1)*a(:,t+1);
    end
    P(:,:,t+1)  = T*P(:,:,t)*transpose(T)-T*P(:,mf,t)*transpose(K(:,:,t)) + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
end
if t<smpl
    PZI_s = PZI;
    K_s = K(:,:,t);
    iF_s = iF(:,:,t);
    P_s = P(:,:,t+1);
    t_steady = t+1;
    P  = cat(3,P(:,:,1:t),repmat(P(:,:,t),[1 1 smpl-t_steady+1]));
    iF = cat(3,iF(:,:,1:t),repmat(inv(P_s(mf,mf)),[1 1 smpl-t_steady+1]));
    L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t_steady+1]));
    K  = cat(3,K(:,:,1:t),repmat(T*P_s(:,mf)*iF_s,[1 1 smpl-t_steady+1]));
end
while t<smpl
    t=t+1;
    v(:,t) = Y(:,t) - a(mf,t) - trend(:,t);
    atilde(:,t) = a(:,t) + PZI*v(:,t);
    a(:,t+1) = T*a(:,t) + K_s*v(:,t);
    aK(1,:,t+1) = a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) = T^(jnk-1)*a(:,t+1);
    end
end
t = smpl+1;
while t>d+1
    t = t-1;
    r(:,t) = Z'*iF(:,:,t)*v(:,t) + L(:,:,t)'*r(:,t+1);
    alphahat(:,t)       = a(:,t) + P(:,:,t)*r(:,t);
    etahat(:,t)         = QRt*r(:,t);
end
if d
    r0 = zeros(mm,d+1); 
    r0(:,d+1) = r(:,d+1);
    r1 = zeros(mm,d+1);
    for t = d:-1:1
        r0(:,t) = Linf(:,:,t)'*r0(:,t+1);
        r1(:,t) = Z'*(iFinf(:,:,t)*v(:,t)-Kstar(:,:,t)'*r0(:,t+1)) + Linf(:,:,t)'*r1(:,t+1);
        alphahat(:,t)   = a(:,t) + Pstar(:,:,t)*r0(:,t) + Pinf(:,:,t)*r1(:,t);
        etahat(:,t)     = QRt*r0(:,t);
    end
end
