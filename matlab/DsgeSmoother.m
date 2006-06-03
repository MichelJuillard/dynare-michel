function [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK] = DsgeSmoother(xparam1,gend,Y)
% stephane.adjemian@cepremap.cnrs.fr [09-07-2004]
%
% Adapted from mj_optmumlik.m
  global bayestopt_ M_ oo_ estim_params_ options_

  alphahat 	= [];
  epsilonhat	= [];
  etahat		= [];
  nobs 		= size(options_.varobs,1);
  smpl        = size(Y,2);

  set_all_parameters(xparam1);

  %------------------------------------------------------------------------------
  % 2. call model setup & reduction program
  %------------------------------------------------------------------------------
  [T,R,SteadyState] = dynare_resolve;
  bayestopt_.mf = bayestopt_.mf2;
  if options_.loglinear == 1
    constant = log(SteadyState(bayestopt_.mfys));
  else
    constant = SteadyState(bayestopt_.mfys);
  end
  trend_coeff = zeros(nobs,1);
  if bayestopt_.with_trend == 1
    trend_coeff = zeros(nobs,1);
    nx1 = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn;
    for i=1:nobs
      trend_coeff(i) = evalin('base',bayestopt_.trend_coeff{i});
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
    Pstar = lyapunov_symm(T,R*Q*transpose(R));
    Pinf	= [];
  elseif options_.lik_init == 2 % Old Diffuse Kalman filter
    Pstar = 10*eye(np);
    Pinf	= [];
  elseif options_.lik_init == 3 % Diffuse Kalman filter
    Pstar = zeros(np,np);
    ivs = bayestopt_.i_T_var_stable;
    Pstar(ivs,ivs) = lyapunov_symm(T(ivs,ivs),R(ivs,:)*Q* ...
				   transpose(R(ivs,:)));
    Pinf  = bayestopt_.Pinf;
    % by M. Ratto
    RR=T(:,find(~ismember([1:np],ivs)));
    i=find(abs(RR)>1.e-10);
    R0=zeros(size(RR));
    R0(i)=sign(RR(i));
    Pinf=R0*R0';
    % by M. Ratto
  end
  % -----------------------------------------------------------------------------
  %  4. Kalman smoother
  % -----------------------------------------------------------------------------
  if estim_params_.nvn
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
	[alphahat,etahat,ahat,aK] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      end
    elseif options_.kalman_algo == 3
      [alphahat,etahat,ahat,aK] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
    end
  end
