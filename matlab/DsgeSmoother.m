function [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T,R,P,PK,d,decomp] = DsgeSmoother(xparam1,gend,Y)
% Estimation of the smoothed variables and innovations. 
% 
% INPUTS 
%   o xparam1       [double]   (p*1) vector of (estimated) parameters. 
%   o gend          [integer]  scalar specifying the number of observations ==> varargin{1}.
%   o data          [double]   (T*n) matrix of data.
%  
% OUTPUTS 
%   o alphahat      [double]  (m*T) matrix, smoothed endogenous variables.
%   o etahat        [double]  (r*T) matrix, smoothed structural shocks (r>n is the umber of shocks).
%   o epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
%   o ahat          [double]  (m*T) matrix, one step ahead filtered (endogenous) variables.
%   o SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable.
%   o trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
%   o aK            [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead filtered (endogenous) variables.
%   o T and R       [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
%    P:             3D array of one-step ahead forecast error variance
%                   matrices
%    PK:            4D array of k-step ahead forecast error variance
%                   matrices (meaningless for periods 1:d)
%    d:             number of periods where filter remains in diffuse part
%                  (should be equal to the order of integration of the model)
%    
% ALGORITHM 
%   Diffuse Kalman filter (Durbin and Koopman)       
%
% SPECIAL REQUIREMENTS
%   None
%  
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.

  global bayestopt_ M_ oo_ estim_params_ options_

  alphahat 	= [];
  etahat	= [];
  epsilonhat	= [];
  ahat          = [];
  SteadyState   = [];
  trend_coeff   = [];
  aK            = [];
  T             = [];
  R             = [];
  P             = [];
  PK            = [];
  d             = [];
  decomp        = [];
  nobs 		= size(options_.varobs,1);
  smpl          = size(Y,2);

  set_all_parameters(xparam1);

  %------------------------------------------------------------------------------
  % 2. call model setup & reduction program
  %------------------------------------------------------------------------------
  [T,R,SteadyState] = dynare_resolve;
  bayestopt_.mf = bayestopt_.mf2;
  if options_.noconstant
      constant = zeros(nobs,1);
  else
      if options_.loglinear == 1
          constant = log(SteadyState(bayestopt_.mfys));
      else
          constant = SteadyState(bayestopt_.mfys);
      end
  end
  trend_coeff = zeros(nobs,1);
  if bayestopt_.with_trend == 1
    trend_coeff = zeros(nobs,1);
    t = options_.trend_coeffs;
    for i=1:length(t)
      if ~isempty(t{i})
	trend_coeff(i) = evalin('base',t{i});
      end
    end
    trend = constant*ones(1,gend)+trend_coeff*(1:gend);
  else
    trend = constant*ones(1,gend);
  end
  start = options_.presample+1;
  np    = size(T,1);
  mf    = bayestopt_.mf;
  % ------------------------------------------------------------------------------
  %  3. Initial condition of the Kalman filter
  % ------------------------------------------------------------------------------
  % 
  %  C'est ici qu'il faut déterminer Pinf et Pstar. Si le modèle est stationnaire,
  %  alors il suffit de poser Pstar comme la solution de l'éuation de Lyapounov et
  %  Pinf=[].
  %
  Q = M_.Sigma_e;
  H = M_.H;
  
  if options_.lik_init == 1		% Kalman filter
    Pstar = lyapunov_symm(T,R*Q*transpose(R),options_.qz_criterium);
    Pinf	= [];
  elseif options_.lik_init == 2 % Old Diffuse Kalman filter
    Pstar = 10*eye(np);
    Pinf	= [];
  elseif options_.lik_init == 3 % Diffuse Kalman filter
    if options_.kalman_algo < 4
      Pstar = zeros(np,np);
      ivs = bayestopt_.restrict_var_list_stationary;
      R1 = R(ivs,:);
      Pstar(ivs,ivs) = lyapunov_symm(T(ivs,ivs),R1*Q*R1',options_.qz_criterium);
      %    Pinf  = bayestopt_.Pinf;
      % by M. Ratto
      RR=T(:,bayestopt_.restrict_var_list_nonstationary);
      i=find(abs(RR)>1.e-10);
      R0=zeros(size(RR));
      R0(i)=sign(RR(i));
      Pinf=R0*R0';
      % by M. Ratto
    else
      [QT,ST] = schur(T);
      e1 = abs(ordeig(ST)) > 2-options_.qz_criterium;
      [QT,ST] = ordschur(QT,ST,e1);
      k = find(abs(ordeig(ST)) > 2-options_.qz_criterium);
      nk = length(k);
      nk1 = nk+1;
      Pinf = zeros(np,np);
      Pinf(1:nk,1:nk) = eye(nk);
      Pstar = zeros(np,np);
      B = QT'*R*Q*R'*QT;
      for i=np:-1:nk+2
	if ST(i,i-1) == 0
	  if i == np
	    c = zeros(np-nk,1);
	  else
	    c = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i,i+1:end)')+...
		ST(i,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i);
	  end
	  q = eye(i-nk)-ST(nk1:i,nk1:i)*ST(i,i);
	  Pstar(nk1:i,i) = q\(B(nk1:i,i)+c);
	  Pstar(i,nk1:i-1) = Pstar(nk1:i-1,i)';
	else
	  if i == np
	    c = zeros(np-nk,1);
	    c1 = zeros(np-nk,1);
	  else
	    c = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i,i+1:end)')+...
		ST(i,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i)+...
		ST(i,i-1)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i-1);
	    c1 = ST(nk1:i,:)*(Pstar(:,i+1:end)*ST(i-1,i+1:end)')+...
		 ST(i-1,i-1)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i-1)+...
		 ST(i-1,i)*ST(nk1:i,i+1:end)*Pstar(i+1:end,i);
	  end
	  q = [eye(i-nk)-ST(nk1:i,nk1:i)*ST(i,i) -ST(nk1:i,nk1:i)*ST(i,i-1);...
	       -ST(nk1:i,nk1:i)*ST(i-1,i) eye(i-nk)-ST(nk1:i,nk1:i)*ST(i-1,i-1)];
	  z =  q\[B(nk1:i,i)+c;B(nk1:i,i-1)+c1];
	  Pstar(nk1:i,i) = z(1:(i-nk));
	  Pstar(nk1:i,i-1) = z(i-nk+1:end);
	  Pstar(i,nk1:i-1) = Pstar(nk1:i-1,i)';
	  Pstar(i-1,nk1:i-2) = Pstar(nk1:i-2,i-1)';
	  i = i - 1;
	end
      end
      if i == nk+2
	c = ST(nk+1,:)*(Pstar(:,nk+2:end)*ST(nk1,nk+2:end)')+ST(nk1,nk1)*ST(nk1,nk+2:end)*Pstar(nk+2:end,nk1);
	Pstar(nk1,nk1)=(B(nk1,nk1)+c)/(1-ST(nk1,nk1)*ST(nk1,nk1));
      end
      
      Z = QT(mf,:);
      R1 = QT'*R;
    end
  end
  % -----------------------------------------------------------------------------
  %  4. Kalman smoother
  % -----------------------------------------------------------------------------
  if any(any(H ~= 0))   % should be replaced by a flag
    if options_.kalman_algo == 1
      [alphahat,epsilonhat,etahat,ahat,aK] = DiffuseKalmanSmootherH1(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      if all(alphahat(:)==0)
	[alphahat,epsilonhat,etahat,ahat,aK] = DiffuseKalmanSmootherH3(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      end
    elseif options_.kalman_algo == 3
      [alphahat,epsilonhat,etahat,ahat,aK] = DiffuseKalmanSmootherH3(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
    end
  else
    if options_.kalman_algo == 1
      [alphahat,etahat,ahat,aK] = DiffuseKalmanSmoother1(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      if all(alphahat(:)==0)
          options_.kalman_algo = 3;
          [alphahat,etahat,ahat,aK] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      end
    elseif options_.kalman_algo == 3
      [alphahat,etahat,ahat,aK] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
    elseif options_.kalman_algo == 4 | options_.kalman_algo == 5
      data1 = Y - trend;
      if options_.kalman_algo == 4
	  [alphahat,etahat,ahat,P,aK,PK,d,decomp] = DiffuseKalmanSmoother1_Z(ST, ...
						  Z,R1,Q,Pinf,Pstar,data1,nobs,np,smpl);
          if all(alphahat(:)==0)
              options_.kalman_algo = 5;
              [alphahat,etahat,ahat,P,aK,PK,d,decomp] = DiffuseKalmanSmoother3_Z(ST, ...
						  Z,R1,Q,Pinf,Pstar, ...
                                                  data1,nobs,np,smpl);
          end
      else
	  [alphahat,etahat,ahat,P,aK,PK,d,decomp] = DiffuseKalmanSmoother3_Z(ST, ...
						  Z,R1,Q,Pinf,Pstar,data1,nobs,np,smpl);
      end
      alphahat = QT*alphahat;
      ahat = QT*ahat;
      if options_.nk > 0
	  nk = options_.nk;
	  for jnk=1:nk
	      aK(jnk,:,:) = QT*squeeze(aK(jnk,:,:));
	      for i=1:size(PK,4)
		  PK(jnk,:,:,i) = QT*squeeze(PK(jnk,:,:,i))*QT';
	      end
	      for i=1:size(decomp,4)
		  decomp(jnk,:,:,i) = QT*squeeze(decomp(jnk,:,:,i));
	      end
	  end
	  for i=1:size(P,4)
	      P(:,:,i) = QT*squeeze(P(:,:,i))*QT';
	  end
      end
    end
  end
