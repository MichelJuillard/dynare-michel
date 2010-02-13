function [alphahat,etahat,a,P,aK,PK,d,decomp] = DiffuseKalmanSmoother3(T,R,Q,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% function [alphahat,etahat,a1, aK] = DiffuseKalmanSmoother3(T,R,Q,Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf)
% Computes the diffuse kalman smoother without measurement error, in the case of a singular var-cov matrix.
% Univariate treatment of multivariate time series.
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
%    a:   matrix of updated variables   (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t})
%
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

% Modified by M. Ratto
% New output argument aK: 1-step to nk-stpe ahed predictions)
% New input argument nk: max order of predictions in aK
% New option options_.diffuse_d where the DKF stops (common with
% diffuselikelihood3)
% New icc variable to count number of iterations for Finf steps
% Pstar % Pinf simmetric
% New termination of DKF iterations based on options_.diffuse_d 
% Li also stored during DKF iterations !!
% some bugs corrected in the DKF part of the smoother (Z matrix and
% alphahat)

global options_

d=0;
decomp=[];
nk = options_.nk;
spinf           = size(Pinf1);
spstar          = size(Pstar1);
v               = zeros(pp,smpl);
a               = zeros(mm,smpl);
a1              = zeros(mm,smpl+1);
aK              = zeros(nk,mm,smpl+nk);
PK              = zeros(nk,mm,mm,smpl+nk);

if isempty(options_.diffuse_d),
    smpl_diff = 1;
else
    smpl_diff=rank(Pinf1);
end

Fstar           = zeros(pp,smpl_diff);
Finf            = zeros(pp,smpl_diff);
Fi              = zeros(pp,smpl_diff);
Ki              = zeros(mm,pp,smpl);
Li              = zeros(mm,mm,pp,smpl);
Linf            = zeros(mm,mm,pp,smpl_diff);
L0              = zeros(mm,mm,pp,smpl_diff);
Kstar           = zeros(mm,pp,smpl_diff);
P               = zeros(mm,mm,smpl+1);
P1                      = P;
Pstar           = zeros(spstar(1),spstar(2),smpl_diff+1); Pstar(:,:,1) = Pstar1;
Pinf            = zeros(spinf(1),spinf(2),smpl_diff+1); Pinf(:,:,1) = Pinf1;
Pstar1          = Pstar;
Pinf1           = Pinf;
crit            = options_.kalman_tol;
crit1       = 1.e-6;
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
icc=0;
newRank   = rank(Pinf(:,:,1),crit1);
while newRank & t < smpl
    t = t+1;
    a(:,t) = a1(:,t);
    Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
    Pinf(:,:,t)=tril(Pinf(:,:,t))+transpose(tril(Pinf(:,:,t),-1));
    Pstar1(:,:,t) = Pstar(:,:,t);
    Pinf1(:,:,t) = Pinf(:,:,t);
    for i=1:pp
        v(i,t)  = Y(i,t)-a(mf(i),t)-trend(i,t);
        Fstar(i,t)      = Pstar(mf(i),mf(i),t);
        Finf(i,t)       = Pinf(mf(i),mf(i),t);
        Kstar(:,i,t)    = Pstar(:,mf(i),t);
        if Finf(i,t) > crit & newRank,  % original MJ: if Finf(i,t) > crit
            icc=icc+1;
            Kinf(:,i,t) = Pinf(:,mf(i),t);
            Linf(:,:,i,t)       = eye(mm) - Kinf(:,i,t)*Z(i,:)/Finf(i,t);
            L0(:,:,i,t)         = (Kinf(:,i,t)*Fstar(i,t)/Finf(i,t) - Kstar(:,i,t))*Z(i,:)/Finf(i,t);
            a(:,t)              = a(:,t) + Kinf(:,i,t)*v(i,t)/Finf(i,t);
            Pstar(:,:,t)        = Pstar(:,:,t) + ...
                Kinf(:,i,t)*transpose(Kinf(:,i,t))*Fstar(i,t)/(Finf(i,t)*Finf(i,t)) - ...
                (Kstar(:,i,t)*transpose(Kinf(:,i,t)) +...
                 Kinf(:,i,t)*transpose(Kstar(:,i,t)))/Finf(i,t);
            Pinf(:,:,t) = Pinf(:,:,t) - Kinf(:,i,t)*transpose(Kinf(:,i,t))/Finf(i,t);
            Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
            Pinf(:,:,t)=tril(Pinf(:,:,t))+transpose(tril(Pinf(:,:,t),-1));
            % new terminiation criteria by M. Ratto
            P0=Pinf(:,:,t);
            %             newRank = any(diag(P0(mf,mf))>crit);
            %             if newRank==0, id = i; end,
            if ~isempty(options_.diffuse_d),  
                newRank = (icc<options_.diffuse_d);  
                %if newRank & any(diag(P0(mf,mf))>crit)==0; 
                if newRank & (any(diag(P0(mf,mf))>crit)==0 & rank(P0,crit1)==0); 
                    disp('WARNING!! Change in OPTIONS_.DIFFUSE_D in univariate DKF')
                    options_.diffuse_d = icc;
                    newRank=0;
                end
            else
                %newRank = any(diag(P0(mf,mf))>crit);                 
                newRank = (any(diag(P0(mf,mf))>crit) | rank(P0,crit1));                 
                if newRank==0, 
                    options_.diffuse_d = icc;
                end                    
            end,
            %             if newRank==0, 
            %                 options_.diffuse_d=i;   %this is buggy
            %             end                    
            % end new terminiation criteria by M. Ratto
        elseif Fstar(i,t) > crit 
            %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
            %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
            %% rank(Pinf)=0. [stéphane,11-03-2004].       
            Li(:,:,i,t)    = eye(mm)-Kstar(:,i,t)*Z(i,:)/Fstar(i,t);  % we need to store Li for DKF smoother
            a(:,t)              = a(:,t) + Kstar(:,i,t)*v(i,t)/Fstar(i,t);
            Pstar(:,:,t)        = Pstar(:,:,t) - Kstar(:,i,t)*transpose(Kstar(:,i,t))/Fstar(i,t);
            Pstar(:,:,t)=tril(Pstar(:,:,t))+transpose(tril(Pstar(:,:,t),-1));
        end
    end
    a1(:,t+1)   = T*a(:,t);
    aK(1,:,t+1) = a1(:,t+1);
    for jnk=2:nk
        aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
    end
    Pstar(:,:,t+1)      = T*Pstar(:,:,t)*transpose(T)+ QQ;
    Pinf(:,:,t+1)       = T*Pinf(:,:,t)*transpose(T);
    P0 = Pinf(:,:,t+1);
    if newRank,
        %newRank = any(diag(P0(mf,mf))>crit);
        newRank   = rank(P0,crit1);
    end
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
        v(i,t)  = Y(i,t) - a(mf(i),t) - trend(i,t);
        Fi(i,t) = P(mf(i),mf(i),t);
        Ki(:,i,t) = P(:,mf(i),t);
        if Fi(i,t) > crit
            Li(:,:,i,t)    = eye(mm)-Ki(:,i,t)*Z(i,:)/Fi(i,t);
            a(:,t) = a(:,t) + Ki(:,i,t)*v(i,t)/Fi(i,t);
            P(:,:,t) = P(:,:,t) - Ki(:,i,t)*transpose(Ki(:,i,t))/Fi(i,t);
            P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
        end
    end
    a1(:,t+1)   = T*a(:,t);
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a1(:,t+1);
    for jnk=1:nk,
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
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
        v(i,t)      = Y(i,t) - a(mf(i),t) - trend(i,t);
        if Fi_s(i) > crit
            a(:,t) = a(:,t) + Ki_s(:,i)*v(i,t)/Fi_s(i);
        end
    end
    a1(:,t+1) = T*a(:,t);
    Pf          = P(:,:,t);
    aK(1,:,t+1) = a1(:,t+1);
    for jnk=1:nk
        Pf = T*Pf*T' + QQ;
        PK(jnk,:,:,t+jnk) = Pf;
        if jnk>1
            aK(jnk,:,t+jnk) = T*dynare_squeeze(aK(jnk-1,:,t+jnk-1));
        end
    end
end
ri=zeros(mm,1);
t = smpl+1;
while t>d+1
    t = t-1;
    for i=pp:-1:1
        if Fi(i,t) > crit
            ri = Z(i,:)'/Fi(i,t)*v(i,t)+Li(:,:,i,t)'*ri;
        end
    end
    r(:,t) = ri;
    alphahat(:,t)       = a1(:,t) + P1(:,:,t)*r(:,t);
    etahat(:,t)         = QRt*r(:,t);
    ri = T'*ri;
end
if d
    r0 = zeros(mm,d);
    r0(:,d) = ri;
    r1 = zeros(mm,d);
    for t = d:-1:1
        for i=pp:-1:1
            if Finf(i,t) > crit & ~(t==d & i>options_.diffuse_d),  % use of options_.diffuse_d to be sure of DKF termination
                                                                   %r1(:,t) = transpose(Z)*v(:,t)/Finf(i,t) + ... BUG HERE in transpose(Z)
                r1(:,t) = Z(i,:)'*v(i,t)/Finf(i,t) + ...
                          L0(:,:,i,t)'*r0(:,t) + Linf(:,:,i,t)'*r1(:,t);
                r0(:,t) = Linf(:,:,i,t)'*r0(:,t);
            elseif Fstar(i,t) > crit % step needed whe Finf == 0
                r0(:,t) = Z(i,:)'/Fstar(i,t)*v(i,t)+Li(:,:,i,t)'*r0(:,t);
            end
        end
        alphahat(:,t)   = a1(:,t) + Pstar1(:,:,t)*r0(:,t) + Pinf1(:,:,t)*r1(:,t);
        r(:,t)  = r0(:,t);
        etahat(:,t)             = QRt*r(:,t);
        if t > 1
            r0(:,t-1) = T'*r0(:,t);
            r1(:,t-1) = T'*r1(:,t);
        end
    end
end

if nargout > 7
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = max(d,1):smpl
        ri_d = zeros(mm,1);
        for i=pp:-1:1
            if Fi(i,t) > crit
                ri_d = Z(i,:)'/Fi(i,t)*v(i,t)+Li(:,:,i,t)'*ri_d;
            end
        end
        
        % calculate eta_tm1t
        eta_tm1t = QRt*ri_d;
        % calculate decomposition
        Ttok = eye(mm,mm); 
        for h = 1:nk
            for j=1:rr
                eta=zeros(rr,1);
                eta(j) = eta_tm1t(j);
                decomp(h,:,j,t+h) = Ttok*P1(:,:,t)*Z'*ZRQinv*Z*R*eta;
            end
            Ttok = T*Ttok;
        end
    end
end

