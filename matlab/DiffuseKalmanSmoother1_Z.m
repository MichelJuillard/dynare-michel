function [alphahat,etahat,atilde,P,aK,PK,d,decomp] = DiffuseKalmanSmoother1_Z(T,Z,R,Q,Pinf1,Pstar1,Y,pp,mm,smpl)

% function [alphahat,etahat,a, aK] = DiffuseKalmanSmoother1(T,Z,R,Q,Pinf1,Pstar1,Y,pp,mm,smpl)
% Computes the diffuse kalman smoother without measurement error, in the case of a non-singular var-cov matrix 
%
% INPUTS
%    T:        mm*mm matrix
%    Z:        pp*mm matrix  
%    R:        mm*rr matrix
%    Q:        rr*rr matrix
%    Pinf1:    mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar1:   mm*mm variance-covariance matrix with stationary variables
%    Y:        pp*1 vector
%    pp:       number of observed variables
%    mm:       number of state variables
%    smpl:     sample size
%             
% OUTPUTS
%    alphahat: smoothed variables (a_{t|T})
%    etahat:   smoothed shocks
%    atilde:   matrix of filtered variables (a_{t|t})
%    aK:       3D array of k step ahead filtered state variables (a_{t+k|t)
%              (meaningless for periods 1:d)
%    P:        3D array of one-step ahead forecast error variance
%              matrices
%    PK:       4D array of k-step ahead forecast error variance
%              matrices (meaningless for periods 1:d)
%    d:        number of periods where filter remains in diffuse part
%              (should be equal to the order of integration of the model)
%    decomp:   decomposition of the effect of shocks on filtered values
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

% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric 


global options_

d = 0;
decomp = [];
nk = options_.nk;
spinf   	= size(Pinf1);
spstar  	= size(Pstar1);
v       	= zeros(pp,smpl);
a       	= zeros(mm,smpl+1);
atilde       	= zeros(mm,smpl);
aK              = zeros(nk,mm,smpl+nk);
PK              = zeros(nk,mm,mm,smpl+nk);
iF      	= zeros(pp,pp,smpl);
Fstar   	= zeros(pp,pp,smpl);
iFinf   	= zeros(pp,pp,smpl);
K       	= zeros(mm,pp,smpl);
L       	= zeros(mm,mm,smpl);
Linf    	= zeros(mm,mm,smpl);
Kstar   	= zeros(mm,pp,smpl);
P       	= zeros(mm,mm,smpl+1);
Pstar   	= zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf    	= zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
crit    	= options_.kalman_tol;
crit1       = 1.e-8;
steady  	= smpl;
rr      	= size(Q,1);
QQ      	= R*Q*transpose(R);
QRt			= Q*transpose(R);
alphahat   	= zeros(mm,smpl);
etahat	   	= zeros(rr,smpl);
r 		   	= zeros(mm,smpl);

t = 0;
while rank(Pinf(:,:,t+1),crit1) & t<smpl
    t = t+1;
    v(:,t)= Y(:,t) - Z*a(:,t);
    F = Z*Pinf(:,:,t)*Z';
    if rcond(F) < crit
    	return		
    end
    iFinf(:,:,t) 	= inv(F);
    PZI                 = Pinf(:,:,t)*Z'*iFinf(:,:,t);
    atilde(:,t)         = a(:,t) + PZI*v(:,t);
    Kinf(:,:,t)	 	= T*PZI;
    a(:,t+1) 	 	= T*atilde(:,t);
    aK(1,:,t+1) 	= a(:,t+1);
    % isn't a meaningless as long as we are in the diffuse part? MJ
    for jnk=2:nk,
        aK(jnk,:,t+jnk) 	 	= T^(jnk-1)*a(:,t+1);
    end
    Linf(:,:,t)  	= T - Kinf(:,:,t)*Z;
    Fstar(:,:,t) 	= Z*Pstar(:,:,t)*Z';
    Kstar(:,:,t) 	= (T*Pstar(:,:,t)*Z'-Kinf(:,:,t)*Fstar(:,:,t))*iFinf(:,:,t);
    Pstar(:,:,t+1)	= T*Pstar(:,:,t)*T'-T*Pstar(:,:,t)*Z'*Kinf(:,:,t)'-Kinf(:,:,t)*F*Kstar(:,:,t)' + QQ;
    Pinf(:,:,t+1)	= T*Pinf(:,:,t)*T'-T*Pinf(:,:,t)*Z'*Kinf(:,:,t)';
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
    v(:,t)      = Y(:,t) - Z*a(:,t);
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    F = Z*P(:,:,t)*Z';
    if rcond(F) < crit
    	return		
    end    
    iF(:,:,t)   = inv(F);
    PZI         = P(:,:,t)*Z'*iF(:,:,t);
    atilde(:,t) = a(:,t) + PZI*v(:,t);
    K(:,:,t)    = T*PZI;
    L(:,:,t)    = T-K(:,:,t)*Z;
    a(:,t+1)    = T*atilde(:,t);
    Pf          = P(:,:,t);
    for jnk=1:nk,
	Pf = T*Pf*T' + QQ;
        aK(jnk,:,t+jnk) = T^jnk*atilde(:,t);
	PK(jnk,:,:,t+jnk) = Pf;
    end
    P(:,:,t+1)  = T*P(:,:,t)*T'-T*P(:,:,t)*Z'*K(:,:,t)' + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
end
if t<smpl
    PZI_s = PZI;
    K_s = K(:,:,t);
    iF_s = iF(:,:,t);
    P_s = P(:,:,t+1);
    P  = cat(3,P(:,:,1:t),repmat(P_s,[1 1 smpl-t]));
    iF = cat(3,iF(:,:,1:t),repmat(iF_s,[1 1 smpl-t]));
    L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t]));
    K  = cat(3,K(:,:,1:t),repmat(T*P_s*Z'*iF_s,[1 1 smpl-t]));
end
while t<smpl
    t=t+1;
    v(:,t) = Y(:,t) - Z*a(:,t);
    atilde(:,t) = a(:,t) + PZI*v(:,t);
    a(:,t+1) = T*atilde(:,t);
    Pf          = P(:,:,t);
    for jnk=1:nk,
	Pf = T*Pf*T' + QQ;
        aK(jnk,:,t+jnk) = T^jnk*atilde(:,t);
	PK(jnk,:,:,t+jnk) = Pf;
    end
end
t = smpl+1;
while t>d+1 & t>2
  t = t-1;
  r(:,t-1) = Z'*iF(:,:,t)*v(:,t) + L(:,:,t)'*r(:,t);
  alphahat(:,t)	= a(:,t) + P(:,:,t)*r(:,t-1);
  etahat(:,t)		= QRt*r(:,t);
end
if d
  r0 = zeros(mm,d); 
  r0(:,d) = r(:,d);
  r1 = zeros(mm,d);
  for t = d:-1:2
    r0(:,t-1) = Linf(:,:,t)'*r0(:,t);
    r1(:,t-1) = Z'*(iFinf(:,:,t)*v(:,t)-Kstar(:,:,t)'*r0(:,t)) + Linf(:,:,t)'*r1(:,t);
    alphahat(:,t)	= a(:,t) + Pstar(:,:,t)*r0(:,t-1) + Pinf(:,:,t)*r1(:,t-1);
    etahat(:,t)		= QRt*r0(:,t);
  end
  r0_0 = Linf(:,:,1)'*r0(:,1);
  r1_0 = Z'*(iFinf(:,:,1)*v(:,1)-Kstar(:,:,1)'*r0(:,1)) + Linf(:,:,1)'*r1(:,1);
  alphahat(:,1)  	= a(:,1) + Pstar(:,:,1)*r0_0 + Pinf(:,:,1)*r1_0;
  etahat(:,1)		= QRt*r0(:,1);
else
  r0 = Z'*iF(:,:,1)*v(:,1) + L(:,:,1)'*r(:,1);
  alphahat(:,1)	= a(:,1) + P(:,:,1)*r0;
  etahat(:,1)	= QRt*r(:,1);
end

if nargout > 7
    decomp = zeros(nk,mm,rr,smpl+nk);
    ZRQinv = inv(Z*QQ*Z');
    for t = d:smpl
        ri_d = Z'*iF(:,:,t)*v(:,t);
        
        % calculate eta_tm1t
        eta_tm1t = QRt*ri_d;
        % calculate decomposition
        Ttok = eye(mm,mm); 
        for h = 1:nk
            for j=1:rr
                eta=zeros(rr,1);
                eta(j) = eta_tm1t(j);
                decomp(h,:,j,t+h) = T^(h-1)*P(:,:,t)*Z'*ZRQinv*Z*R*eta;
            end
        end
    end
end
