function [LIK, lik] = DiffuseLikelihoodH3(T,R,Q,H,Pinf,Pstar,Y,trend,start)

% function [LIK, lik] = DiffuseLikelihoodH3(T,R,Q,H,Pinf,Pstar,Y,trend,start)
% Computes the diffuse likelihood with measurement error, in the case of
% a singular var-cov matrix.
% Univariate treatment of multivariate time series.
%
% INPUTS
%    T:      mm*mm matrix
%    R:      mm*rr matrix
%    Q:      rr*rr matrix
%    H:      pp*pp matrix
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
% part of DYNARE, copyright Dynare Team (2005-2008)
% Gnu Public License.



% M. Ratto added lik in output [October 2005]
% changes by M. Ratto
% introduced new global variable id_ for termination of DKF
% introduced a persistent fmax, in order to keep track the max order of
% magnitude of the 'zero' values in Pinf at DKF termination
% new icc counter for Finf steps in DKF
% new termination for DKF
% likelihood terms for Fstar must be cumulated in DKF also when Pinf is non
% zero. this bug is fixed.

global bayestopt_ options_
  
mf = bayestopt_.mf;
pp     = size(Y,1);
mm     = size(T,1);
smpl   = size(Y,2);
a      = zeros(mm,1);
QQ     = R*Q*transpose(R);
t      = 0;
lik    = zeros(smpl+1,1);
lik(smpl+1) = smpl*pp*log(2*pi); %% the constant of minus two times the log-likelihood 
notsteady 	= 1;
crit      	= options_.kalman_tol;
crit1      	= 1.e-6;
newRank	  	= rank(Pinf,crit1);
icc = 0;
while newRank & t < smpl %% Matrix Finf is assumed to be zero
  t = t+1;
  for i=1:pp
    v(i) 	= Y(i,t)-a(mf(i))-trend(i,t);
    Fstar 	= Pstar(mf(i),mf(i))+H(i,i);
    Finf	= Pinf(mf(i),mf(i));
    Kstar 	= Pstar(:,mf(i));
    if Finf > crit & newRank
      icc = icc + 1;
      Kinf	= Pinf(:,mf(i));
      a		= a + Kinf*v(i)/Finf;
      Pstar	= Pstar + Kinf*transpose(Kinf)*Fstar/(Finf*Finf) - ...
	  (Kstar*transpose(Kinf)+Kinf*transpose(Kstar))/Finf;
      Pinf	= Pinf - Kinf*transpose(Kinf)/Finf;
      lik(t) 	= lik(t) + log(Finf);
      % start new termination criterion for DKF
      if ~isempty(options_.diffuse_d),  
	newRank = (icc<options_.diffuse_d);  
	%if newRank & any(diag(Pinf(mf,mf))>crit)==0; %  M. Ratto this line is BUGGY
	if newRank & (any(diag(Pinf(mf,mf))>crit)==0 & rank(Pinf,crit1)==0); 
	  options_.diffuse_d = icc;
	  newRank=0;
	  disp('WARNING: Change in OPTIONS_.DIFFUSE_D in univariate DKF')
	  disp(['new OPTIONS_.DIFFUSE_D = ',int2str(icc)])
	  disp('You may have to reset the optimisation')
	end
      else
	%newRank = any(diag(Pinf(mf,mf))>crit);     % M. Ratto this line is BUGGY 
	newRank = (any(diag(Pinf(mf,mf))>crit) | rank(Pinf,crit1));                 
	if newRank==0, 
	  P0=	T*Pinf*transpose(T);
	  %newRank = any(diag(P0(mf,mf))>crit);   % M. Ratto this line is BUGGY
	  newRank = (any(diag(P0(mf,mf))>crit) | rank(P0,crit1));   % M. Ratto 10 Oct 2005
	  if newRank==0, 
	    options_.diffuse_d = icc;
	  end
	end                    
      end,
      % end new termination and checks for DKF and fmax
    elseif Finf > crit 
      %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
      %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
      %% rank(Pinf)=0. [stéphane,11-03-2004].	  
      %if rank(Pinf) == 0
      % the likelihood terms should alwasy be cumulated, not only
      % when Pinf=0, otherwise the lik would depend on the ordering
      % of observed variables
      lik(t)	= lik(t) + log(Fstar) + v(i)*v(i)/Fstar;
      %end
      a 	= a + Kstar*v(i)/Fstar;
      Pstar	= Pstar - Kstar*transpose(Kstar)/Fstar;					
    else
      % disp(['zero F term in DKF for observed ',int2str(i),' ',num2str(Fi)])
    end
  end
  if newRank
    oldRank = rank(Pinf,crit1);
  else
    oldRank = 0;
  end
  a 		= T*a;
  Pstar 	= T*Pstar*transpose(T)+QQ;
  Pinf	= T*Pinf*transpose(T);
  if newRank
    newRank = rank(Pinf,crit1);
  end
  if oldRank ~= newRank
    disp('DiffuseLiklihoodH3 :: T does influence the rank of Pinf!')	
  end		 		
end
if t == smpl                                                           
  error(['There isn''t enough information to estimate the initial' ... 
	 ' conditions of the nonstationary variables']);                   
end                                                                    
while notsteady & t < smpl
  t = t+1;
  for i=1:pp
    v(i) = Y(i,t) - a(mf(i)) - trend(i,t);
    Fi   = Pstar(mf(i),mf(i))+H(i,i);
    if Fi > crit
      Ki	= Pstar(:,mf(i));
      a		= a + Ki*v(i)/Fi;
      Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
      lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
    end
  end	
  oldP 		= Pstar;
  a 		= T*a;
  Pstar 	= T*Pstar*transpose(T) + QQ;
  notsteady 	= ~(max(max(abs(Pstar-oldP)))<crit);
end
while t < smpl
  t = t+1;
  for i=1:pp
    v(i) = Y(i,t) - a(mf(i)) - trend(i,t);
    Fi   = Pstar(mf(i),mf(i))+H(i,i);
    if Fi > crit
      Ki 		= Pstar(:,mf(i));
      a 		= a + Ki*v(i)/Fi;
      Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
      lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
    end
  end	
  a = T*a;
end
LIK = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);

