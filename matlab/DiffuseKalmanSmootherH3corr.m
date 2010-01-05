function [alphahat,epsilonhat,etahat,a,aK] = DiffuseKalmanSmootherH3corr(T,R,Q,H,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)

% function [alphahat,epsilonhat,etahat,a1] = DiffuseKalmanSmootherH3corr(T,R,Q,H,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% Computes the diffuse kalman smoother with measurement error, in the case of a singular var-cov matrix.
% Univariate treatment of multivariate time series.
%
% INPUTS
%    T:         mm*mm matrix
%    R:         mm*rr matrix
%    Q:         rr*rr matrix
%    Pinf1:     mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:    mm*mm variance-covariance matrix with stationary variables
%    Y:         pp*1 vector
%    trend
%    pp:        number of observed variables
%    mm:        number of state variables
%    smpl:      sample size
%    mf:        observed variables index in the state vector
%             
% OUTPUTS
%    alphahat:  smoothed state variables (a_{t|T})
%    epsilonhat:smoothed measurement error
%    etahat:    smoothed shocks
%    a:         matrix of updated variables (a_{t|t})
%    aK:        matrix of one step ahead filtered state variables (a_{t+k|t})

% SPECIAL REQUIREMENTS
%   See "Fast Filtering and Smoothing for Multivariate State Space
%   Models", S.J. Koopman and J. Durbin (2000, in Journal of Time Series 
%   Analysis, vol. 21(3), pp. 281-296).  

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

global options_;

nk = options_.nk;
rr = size(Q,1);
T  = cat(1,cat(2,T,zeros(mm,pp)),zeros(pp,mm+pp));
R  = cat(1,cat(2,R,zeros(mm,pp)),cat(2,zeros(pp,rr),eye(pp)));
Q  = cat(1,cat(2,Q,zeros(rr,pp)),cat(2,zeros(pp,rr),H));
if size(Pinf1,1) % Otherwise Pinf = 0 (no unit root) 
    Pinf1 = cat(1,cat(2,Pinf1,zeros(mm,pp)),zeros(pp,mm+pp));
end
Pstar1   = cat(1,cat(2,Pstar1,zeros(mm,pp)),cat(2,zeros(pp,mm),H));
spinf    = size(Pinf1);
spstar   = size(Pstar1);
Pstar    = zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf     = zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
Pstar1   = Pstar;
Pinf1    = Pinf;
v        = zeros(pp,smpl);
a        = zeros(mm+pp,smpl);
a1       = zeros(mm+pp,smpl+1);
aK       = zeros(nk,mm,smpl+nk);
Fstar    = zeros(pp,smpl);
Finf     = zeros(pp,smpl);
Fi       = zeros(pp,smpl);
Ki       = zeros(mm+pp,pp,smpl);
Li       = zeros(mm+pp,mm+pp,pp,smpl);
Linf     = zeros(mm+pp,mm+pp,pp,smpl);
L0       = zeros(mm+pp,mm+pp,pp,smpl);
Kstar    = zeros(mm+pp,pp,smpl);
Kinf     = zeros(mm+pp,pp,smpl);
P        = zeros(mm+pp,mm+pp,smpl+1);
P1               = zeros(mm+pp,mm+pp,smpl+1);
crit     = options_.kalman_tol;
steady   = smpl;
QQ       = R*Q*transpose(R);
QRt              = Q*transpose(R);
alphahat        = zeros(mm+pp,smpl);
etahat          = zeros(rr,smpl);
epsilonhat      = zeros(pp,smpl);
r                       = zeros(mm+pp,smpl+1);
Z = zeros(pp,mm+pp);
for i=1:pp;
    Z(i,mf(i)) = 1;
    Z(i,mm+i)  = 1;
end
%% [1] Kalman filter...
t = 0;
newRank   = rank(Pinf(:,:,1),crit);
while newRank & t < smpl
    t = t+1;
    a(:,t) = a1(:,t);
    Pstar1(:,:,t) = Pstar(:,:,t);
    Pinf1(:,:,t) = Pinf(:,:,t);
    for i=1:pp
        v(i,t)  = Y(i,t)-a(mf(i),t)-a(mm+i,t)-trend(i,t);
        Fstar(i,t)       = Pstar(mf(i),mf(i),t)+Pstar(mm+i,mm+i,t);
        Finf(i,t)        = Pinf(mf(i),mf(i),t);
        Kstar(:,i,t) = Pstar(:,mf(i),t)+Pstar(:,mm+i,t);
        if Finf(i,t) > crit
            Kinf(:,i,t) = Pinf(:,mf(i),t);
            Linf(:,:,i,t)       = eye(mm+pp) - Kinf(:,i,t)*Z(i,:)/Finf(i,t);
            L0(:,:,i,t)         = (Kinf(:,i,t)*Fstar(i,t)/Finf(i,t) - Kstar(:,i,t))*Z(i,:)/Finf(i,t);
            a(:,t)              = a(:,t) + Kinf(:,i,t)*v(i,t)/Finf(i,t);
            Pstar(:,:,t)        = Pstar(:,:,t) + ...
                Kinf(:,i,t)*transpose(Kinf(:,i,t))*Fstar(i,t)/(Finf(i,t)*Finf(i,t)) - ...
                (Kstar(:,i,t)*transpose(Kinf(:,i,t)) +...
                 Kinf(:,i,t)*transpose(Kstar(:,i,t)))/Finf(i,t);
            Pinf(:,:,t) = Pinf(:,:,t) - Kinf(:,i,t)*transpose(Kinf(:,i,t))/Finf(i,t);
        else %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
            %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
            %% rank(Pinf)=0. [stéphane,11-03-2004].       
            a(:,t)              = a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);
            Pstar(:,:,t)        = Pstar(:,:,t) - Kstar(:,i,t)*transpose(Kstar(:,i,t))/Fstar(i,t);
        end
    end
    a1(:,t+1)           = T*a(:,t);
    for jnk=1:nk,
        aK(jnk,:,t+jnk)                 = T^jnk*a(:,t);
    end
    Pstar(:,:,t+1)      = T*Pstar(:,:,t)*transpose(T)+ QQ;
    Pinf(:,:,t+1)       = T*Pinf(:,:,t)*transpose(T);
    P0=Pinf(:,:,t+1);
    newRank = ~all(abs(P0(:))<crit);
end
d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
Linf  = Linf(:,:,:,1:d);
L0  = L0(:,:,:,1:d);
Fstar = Fstar(:,1:d);
Finf = Finf(:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
Pstar1 = Pstar1(:,:,1:d);
Pinf1  = Pinf1(:,:,1:d);
notsteady = 1;
while notsteady & t<smpl
    t = t+1;
    a(:,t) = a1(:,t);
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    P1(:,:,t) = P(:,:,t);
    for i=1:pp
        v(i,t)    = Y(i,t) - a(mf(i),t) - a(mm+i,t) - trend(i,t);
        Fi(i,t)   = P(mf(i),mf(i),t)+P(mm+i,mm+i,t);
        Ki(:,i,t) = P(:,mf(i),t)+P(:,mm+i,t);
        if Fi(i,t) > crit
            Li(:,:,i,t)    = eye(mm+pp)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
            a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);
            P(:,:,t) = P(:,:,t) - Ki(:,i,t)*transpose(Ki(:,i,t))/Fi(i,t);
            P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
        end
    end
    a1(:,t+1) = T*a(:,t);
    for jnk=1:nk,
        aK(jnk,:,t+jnk)                 = T^jnk*a(:,t);
    end
    P(:,:,t+1) = T*P(:,:,t)*transpose(T) + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
end
P_s=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
Fi_s = Fi(:,t);
Ki_s = Ki(:,:,t);
L_s  =Li(:,:,:,t);
if t<smpl
    t_steady = t+1;
    P  = cat(3,P(:,:,1:t),repmat(P(:,:,t),[1 1 smpl-t_steady+1]));
    Fi = cat(2,Fi(:,1:t),repmat(Fi_s,[1 1 smpl-t_steady+1]));
    Li  = cat(4,Li(:,:,:,1:t),repmat(L_s,[1 1 smpl-t_steady+1]));
    Ki  = cat(3,Ki(:,:,1:t),repmat(Ki_s,[1 1 smpl-t_steady+1]));
end
while t<smpl
    t=t+1;
    a(:,t) = a1(:,t);
    for i=1:pp
        v(i,t)      = Y(i,t) - a(mf(i),t) - a(mm+i,t) - trend(i,t);
        if Fi_s(i) > crit
            a(:,t) = a(:,t) + Ki_s(:,i)*v(i,t)/Fi_s(i);
        end
    end
    a1(:,t+1) = T*a(:,t);
    for jnk=1:nk,
        aK(jnk,:,t+jnk)                 = T^jnk*a(:,t);
    end
end

%% [2] Kalman smoother...
ri=zeros(mm,1);
t = smpl+1;
while t>d+1 
    t = t-1;
    for i=pp:-1:1
        if Fi(i,t) > crit
            ri=Z(i,:)'/Fi(i,t)*v(i,t)+Li(:,:,i,t)'*ri;
        end
    end
    r(:,t) = ri(:,t);
    alphahat(:,t) = a1(:,t) + P1(:,:,t)*r(:,t);
    tmp   = QRt*r(:,t);
    etahat(:,t)         = tmp(1:rr);
    epsilonhat(:,t)     = tmp(rr+(1:pp));
    ri = T'*ri;
end
if d
    r0 = zeros(mm+pp,d); 
    r0(:,d) = ri;
    r1 = zeros(mm+pp,d);
    for t = d:-1:1
        for i=pp:-1:1
            if Finf(i,t) > crit
                r1(:,t) = transpose(Z)*v(:,t)/Finf(i,t) + ...
                          L0(:,:,i,t)'*r0(:,t) + Linf(:,:,i,t)'*r1(:,t);
                r0(:,t) = transpose(Linf(:,:,i,t))*r0(:,t);
            end
        end
        alphahat(:,t) = a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t);
        r(:,t-1)      = r0(:,t);
        tmp              = QRt*r(:,t);
        etahat(:,t)   = tmp(1:rr);
        epsilonhat(:,t)    = tmp(rr+(1:pp));
        if t > 1
            r0(:,t-1) = T'*r0(:,t);
            r1(:,t-1) = T'*r1(:,t);
        end
    end
end
alphahat = alphahat(1:mm,:);