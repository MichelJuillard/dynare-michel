function [LIK, lik] = DiffuseLikelihood1(T,R,Q,Pinf,Pstar,Y,trend,start)

% function [LIK, lik] = DiffuseLikelihood1(T,R,Q,Pinf,Pstar,Y,trend,start)
% Computes the diffuse likelihood without measurement error, in the case of a non-singular var-cov matrix 
%
% INPUTS
%    T:      mm*mm matrix
%    R:      mm*rr matrix
%    Q:      rr*rr matrix
%    Pinf:   mm*mm diagonal matrix with with q ones and m-q zeros
%    Pstar:  mm*mm variance-covariance matrix with stationary variables
%    Y:      pp*1 vector
%    trend
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
%  
% part of DYNARE, copyright Dynare Team (2004-2008)
% Gnu Public License.



% M. Ratto added lik in output

  global bayestopt_ options_
  
  mf = bayestopt_.mf;
  smpl = size(Y,2);
  mm   = size(T,2);
  pp   = size(Y,1);
  a    = zeros(mm,1);
  dF   = 1;
  QQ   = R*Q*transpose(R);
  t    = 0;
  lik  = zeros(smpl+1,1);
  LIK  = Inf;
  lik(smpl+1) = smpl*pp*log(2*pi);
  notsteady   = 1;
  crit        = options_.kalman_tol;
  reste       = 0;
  while rank(Pinf,crit) & t < smpl
    t     = t+1;
    v  	  = Y(:,t)-a(mf)-trend(:,t);
    Finf  = Pinf(mf,mf);
    if rcond(Finf) < crit 
      if ~all(abs(Finf(:)) < crit)
	return
      else
	iFstar	= inv(Pstar(mf,mf));
	dFstar	= det(Pstar(mf,mf));
	Kstar	= Pstar(:,mf)*iFstar;
	lik(t)	= log(dFstar) + transpose(v)*iFstar*v;
	Pinf	= T*Pinf*transpose(T);
	Pstar	= T*(Pstar-Pstar(:,mf)*transpose(Kstar))*transpose(T)+QQ;
	a	= T*(a+Kstar*v);
      end
    else
      lik(t)	= log(det(Finf));
      iFinf	= inv(Finf);
      Kinf	= Pinf(:,mf)*iFinf;			%%	premultiplication by the transition matrix T is removed (stephane) 
      Fstar	= Pstar(mf,mf);
      Kstar	= (Pstar(:,mf)-Kinf*Fstar)*iFinf; 	%%	premultiplication by the transition matrix T is removed (stephane)
      Pstar	= T*(Pstar-Pstar(:,mf)*transpose(Kinf)-Pinf(:,mf)*transpose(Kstar))*transpose(T)+QQ;
      Pinf	= T*(Pinf-Pinf(:,mf)*transpose(Kinf))*transpose(T);
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
    v  	  = Y(:,t)-a(mf)-trend(:,t);
    F  = Pstar(mf,mf);
    oldPstar  = Pstar;
    dF = det(F);
    if rcond(F) < crit 
      if ~all(abs(F(:))<crit)
	return
      else
	a         = T*a;
	Pstar     = T*Pstar*transpose(T)+QQ;
      end
    else
      F_singular = 0;
      iF        = inv(F);
      lik(t)    = log(dF)+transpose(v)*iF*v;
      K         = Pstar(:,mf)*iF; %% premultiplication by the transition matrix T is removed (stephane)
      a         = T*(a+K*v);		%% --> factorization of the transition matrix...
      Pstar     = T*(Pstar-K*Pstar(mf,:))*transpose(T)+QQ;	%% ... idem (stephane)
    end
    notsteady = ~(max(max(abs(Pstar-oldPstar)))<crit);
  end
  if F_singular == 1
    error(['The variance of the forecast error remains singular until the' ...
	  'end of the sample'])
  end
  reste = smpl-t;
  while t < smpl
    t = t+1;
    v = Y(:,t)-a(mf)-trend(:,t);
    a = T*(a+K*v);
    lik(t) = transpose(v)*iF*v;
  end
  lik(t) = lik(t) + reste*log(dF);


  LIK    = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);% Minus the log-likelihood.
