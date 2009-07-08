function [LIK, lik] = DiffuseLikelihood1_Z(T,Z,R,Q,Pinf,Pstar,Y,start)

% function [LIK, lik] = DiffuseLikelihood1_Z(T,Z,R,Q,Pinf,Pstar,Y,start)
% Computes the diffuse likelihood without measurement error, in the case of a non-singular var-cov matrix 
%
% INPUTS
%    T:      mm*mm matrix
%    Z:      pp,mm matrix  
%    R:      mm*rr matrix
%    Q:      rr*rr matrix
%    Pinf:   mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar:  mm*mm variance-covariance matrix with stationary variables
%    Y:      pp*1 vector
%    start:  likelihood evaluation at 'start'
%             
% OUTPUTS
%    LIK:    likelihood
%    lik:    density vector in each period
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

% M. Ratto added lik in output

  global bayestopt_ options_
  
  smpl = size(Y,2);
  mm   = size(T,2);
  pp   = size(Y,1);
  a    = zeros(mm,1);
  dF   = 1;
  QQ   = R*Q*transpose(R);
  t    = 0;
  lik  = zeros(smpl,1);
  LIK  = Inf;
  notsteady   = 1;
  crit        = options_.kalman_tol;
  while rank(Pinf,crit) & t < smpl
    t     = t+1;
    v  	  = Y(:,t)-Z*a;
    Finf  = Z*Pinf*Z';
    if rcond(Finf) < crit 
      if ~all(abs(Finf(:)) < crit)
	return
      else
	Fstar   = Z*Pstar*Z';
	iFstar	= inv(Fstar);
	dFstar	= det(Fstar);
	Kstar	= Pstar*Z'*iFstar;
	lik(t)	= log(dFstar) + v'*iFstar*v;
	Pinf	= T*Pinf*transpose(T);
	Pstar	= T*(Pstar-Pstar*Z'*Kstar')*T'+QQ;
	a	= T*(a+Kstar*v);
      end
    else
      lik(t)	= log(det(Finf));
      iFinf	= inv(Finf);
      Kinf	= Pinf*Z'*iFinf;		
      Fstar	= Z*Pstar*Z';
      Kstar	= (Pstar*Z'-Kinf*Fstar)*iFinf; 	
      Pstar	= T*(Pstar-Pstar*Z'*Kinf'-Pinf*Z'*Kstar')*T'+QQ;
      Pinf	= T*(Pinf-Pinf*Z'*Kinf')*T';
      a		= T*(a+Kinf*v);					
    end  
  end
  if t == smpl                                                           
    error(['There isn''t enough information to estimate the initial' ... 
	   ' conditions of the nonstationary variables']);                   
  end                                                                    
  F_singular = 1;
  while notsteady & t < smpl
    t  = t+1;
    v  	  = Y(:,t)-Z*a;
    F  = Z*Pstar*Z';
    oldPstar  = Pstar;
    dF = det(F);
    if rcond(F) < crit 
      if ~all(abs(F(:))<crit)
	return
      else
	a         = T*a;
	Pstar     = T*Pstar*T'+QQ;
      end
    else
      F_singular = 0;
      iF        = inv(F);
      lik(t)    = log(dF)+v'*iF*v;
      K         = Pstar*Z'*iF;
      a         = T*(a+K*v);	
      Pstar     = T*(Pstar-K*Z*Pstar)*T'+QQ;
    end
    notsteady = ~(max(max(abs(Pstar-oldPstar)))<crit);
  end
  if F_singular == 1
    error(['The variance of the forecast error remains singular until the' ...
	  'end of the sample'])
  end
  if t < smpl
    t0 = t+1;
    while t < smpl
      t = t+1;
      v = Y(:,t)-Z*a;
      a = T*(a+K*v);
      lik(t) = v'*iF*v;
    end
    lik(t0:smpl) = lik(t0:smpl) + log(dF);
  end
  % adding log-likelihhod constants
  lik = (lik + pp*log(2*pi))/2;

  LIK = sum(lik(start:end)); % Minus the log-likelihood.