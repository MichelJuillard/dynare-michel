function [fval,llik,cost_flag,ys,trend_coeff,info] = dsge_likelihood_hh(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% function [fval,llik,cost_flag,ys,trend_coeff,info] = DsgeLikelihood_hh(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% Evaluates the posterior kernel of a dsge model.
%
% INPUTS
%   xparam1                        [double]   vector of model parameters.
%   gend                           [integer]  scalar specifying the number of observations.
%   data                           [double]   matrix of data
%   data_index                     [cell]     cell of column vectors
%   number_of_observations         [integer]
%   no_more_missing_observations   [integer]
% OUTPUTS
%   fval        :     value of the posterior kernel at xparam1.
%   llik        :     probabilities at each time point
%   cost_flag   :     zero if the function returns a penalty, one otherwise.
%   ys          :     steady state of original endogenous variables
%   trend_coeff :
%   info        :     vector of informations about the penalty:
%                     41: one (many) parameter(s) do(es) not satisfied the lower bound
%                     42: one (many) parameter(s) do(es) not satisfied the upper bound
%
% SPECIAL REQUIREMENTS
%

% Copyright (C) 2004-2011 Dynare Team
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


% Declaration of the penalty as a persistent variable.

% Persistent variable 'penalty' is used to compute an endogenous penalty to
% the value 'fval' when various conditions are encountered. These conditions
% set also 'exit_flag' equal to 0 instead of 1.  It is only when
% dsge_likelihood_hh() is called by an newrat() called by
% dynare_estimation_1() that 'exit_flag' is ignored and penalized 'fval' is
% actually used.  
% In that case, 'penalty' is properly initialized, at the very end of the
% present function, by a call to dsge_likelihood_hh() made in
% initial_estimation_checks(). If a condition triggers exit_flag ==
% 0, initial_estimation_checks() triggers an error.
% In summary, an initial call to the present function, without triggering
% any condition, guarantees that 'penalty' is properly initialized when needed.
persistent penalty

% Initialization of the returned variables and others...
fval        = [];
ys          = [];
trend_coeff = [];
cost_flag   = 1;
llik=NaN;
info        = 0;
singularity_flag = 0;

if DynareOptions.block == 1
    error('DsgeLikelihood_hh:: This routine (called if mode_compute==5) is not compatible with the block option!')
end

%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------

% Return, with endogenous penalty, if some parameters are smaller than the lower bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparam1<BayesInfo.lb)
    k = find(xparam1<BayesInfo.lb);
    fval = penalty+sum((BayesInfo.lb(k)-xparam1(k)).^2);
    exit_flag = 0;
    info = 41;
    return
end

% Return, with endogenous penalty, if some parameters are greater than the upper bound of the prior domain.
if ~isequal(DynareOptions.mode_compute,1) && any(xparam1>BayesInfo.ub)
    k = find(xparam1>BayesInfo.ub);
    fval = penalty+sum((xparam1(k)-BayesInfo.ub(k)).^2);
    exit_flag = 0;
    info = 42;
    return
end

% Get the diagonal elements of the covariance matrices for the structural innovations (Q) and the measurement error (H).
Q = Model.Sigma_e;
H = Model.H;
for i=1:EstimatedParameters.nvx
    k =EstimatedParameters.var_exo(i,1);
    Q(k,k) = xparam1(i)*xparam1(i);
end
offset = EstimatedParameters.nvx;
if EstimatedParameters.nvn
    for i=1:EstimatedParameters.nvn
        k = EstimatedParameters.var_endo(i,1);
        H(k,k) = xparam1(i+offset)*xparam1(i+offset);
    end
    offset = offset+EstimatedParameters.nvn;
else
    H = zeros(DynareDataset.info.nvobs);
end

% Get the off-diagonal elements of the covariance matrix for the structural innovations. Test if Q is positive definite.
if EstimatedParameters.ncx
    for i=1:EstimatedParameters.ncx
        k1 =EstimatedParameters.corrx(i,1);
        k2 =EstimatedParameters.corrx(i,2);
        Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
        Q(k2,k1) = Q(k1,k2);
    end
    % Try to compute the cholesky decomposition of Q (possible iff Q is positive definite)
    [CholQ,testQ] = chol(Q);
    if testQ
        % The variance-covariance matrix of the structural innovations is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        a = diag(eig(Q));
        k = find(a < 0);
        if k > 0
            fval = penalty+sum(-a(k));
            exit_flag = 0;
            info = 43;
            return
        end
    end
    offset = offset+EstimatedParameters.ncx;
end

% Get the off-diagonal elements of the covariance matrix for the measurement errors. Test if H is positive definite.
if EstimatedParameters.ncn
    for i=1:EstimatedParameters.ncn
        k1 = DynareOptions.lgyidx2varobs(EstimatedParameters.corrn(i,1));
        k2 = DynareOptions.lgyidx2varobs(EstimatedParameters.corrn(i,2));
        H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
        H(k2,k1) = H(k1,k2);
    end
    % Try to compute the cholesky decomposition of H (possible iff H is positive definite)
    [CholH,testH] = chol(H);
    if testH
        % The variance-covariance matrix of the measurement errors is not definite positive. We have to compute the eigenvalues of this matrix in order to build the endogenous penalty.
        a = diag(eig(H));
        k = find(a < 0);
        if k > 0
            fval = penalty+sum(-a(k));
            exit_flag = 0;
            info = 44;
            return
        end
    end
    offset = offset+EstimatedParameters.ncn;
end

% Update estimated structural parameters in Mode.params.
if EstimatedParameters.np > 0
    Model.params(EstimatedParameters.param_vals(:,1)) = xparam1(offset+1:end);
end

% Update Model.Sigma_e and Model.H.
Model.Sigma_e = Q;
Model.H = H;

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

% Linearize the model around the deterministic sdteadystate and extract the matrices of the state equation (T and R).
[T,R,SteadyState,info,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults,'restrict');

% Return, with endogenous penalty when possible, if dynare_resolve issues an error code (defined in resol).
if info(1) == 1 || info(1) == 2 || info(1) == 5 || info(1) == 22 || info(1) == 24
    fval = penalty+1;
    info = info(1);
    exit_flag = 0;
    return
elseif info(1) == 3 || info(1) == 4 || info(1)==6 ||info(1) == 19 || info(1) == 20 || info(1) == 21  || info(1) == 23
    fval = penalty+info(2);
    info = info(1);
    exit_flag = 0;
    return
end

% Define a vector of indices for the observed variables. Is this really usefull?...
BayesInfo.mf = BayesInfo.mf1;

% Define the constant vector of the measurement equation.
if DynareOptions.noconstant
    constant = zeros(DynareDataset.info.nvobs,1);
else
    if DynareOptions.loglinear
        constant = log(SteadyState(BayesInfo.mfys));
    else
        constant = SteadyState(BayesInfo.mfys);
    end
end

% Define the deterministic linear trend of the measurement equation.
if BayesInfo.with_trend
    trend_coeff = zeros(DynareDataset.info.nvobs,1);
    t = DynareOptions.trend_coeffs;
    for i=1:length(t)
        if ~isempty(t{i})
            trend_coeff(i) = evalin('base',t{i});
        end
    end
    trend = repmat(constant,1,DynareDataset.info.ntobs)+trend_coeff*[1:DynareDataset.info.ntobs];
else
    trend = repmat(constant,1,DynareDataset.info.ntobs);
end

% Get needed informations for kalman filter routines.
start = DynareOptions.presample+1;
Z = BayesInfo.mf; % old mf
no_missing_data_flag = ~DynareDataset.missing.state;
mm = length(T); % old np
pp = DynareDataset.info.nvobs;
rr = length(Q);
kalman_tol = DynareOptions.kalman_tol;
riccati_tol = DynareOptions.riccati_tol;
Y   = DynareDataset.data-trend;

%------------------------------------------------------------------------------
% 3. Initial condition of the Kalman filter
%------------------------------------------------------------------------------
kalman_algo = DynareOptions.kalman_algo;

% resetting measurement errors covariance matrix for univariate filters
if (kalman_algo == 2) || (kalman_algo == 4)
    if isequal(H,0)
        H = zeros(nobs,1);
        mmm = mm;
    else
        if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
            H = diag(H);
            mmm = mm; 
        else
            Z = [Z, eye(pp)];
            T = blkdiag(T,zeros(pp));
            Q = blkdiag(Q,H);
            R = blkdiag(R,eye(pp));
            Pstar = blkdiag(Pstar,H);
            Pinf  = blckdiag(Pinf,zeros(pp));
            H = zeros(nobs,1);
            mmm   = mm+pp;
        end
    end
end


diffuse_periods = 0;
correlated_errors_have_been_checked = 0;
singular_diffuse_filter = 0;
switch DynareOptions.lik_init
  case 1% Standard initialization with the steady state of the state equation.
    if kalman_algo~=2
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
        kalman_algo = 1;
    end
    Pstar = lyapunov_symm(T,R*Q*R',DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
    Pinf  = [];
    a     = zeros(mm,1);
    Zflag = 0;
  case 2% Initialization with large numbers on the diagonal of the covariance matrix if the states (for non stationary models).
    if kalman_algo ~= 2
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
        kalman_algo = 1;
    end
    Pstar = DynareOptions.Harvey_scale_factor*eye(mm);
    Pinf  = [];
    a     = zeros(mm,1);
    Zflag = 0;
  case 3% Diffuse Kalman filter (Durbin and Koopman)
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
    if kalman_algo == 0
        kalman_algo = 3;
    elseif ~((kalman_algo == 3) || (kalman_algo == 4)) 
        error(['diffuse filter: options_.kalman_algo can only be equal ' ...
               'to 0 (default), 3 or 4'])
    end

    [Z,T,R,QT,Pstar,Pinf] = schur_statespace_transformation(Z,T,R,Q,DynareOptions.qz_criterium);
    Zflag = 1;
    % Run diffuse kalman filter on first periods.
    if (kalman_algo==3)
        % Multivariate Diffuse Kalman Filter
        if no_missing_data_flag
            [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, 1, size(Y,2), ...
                                                  zeros(mm,1), Pinf, Pstar, ...
                                                  kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                  T,R,Q,H,Z,mm,pp,rr);
        else
            [dLIK,dlik,a,Pstar] = missing_observations_kalman_filter_d(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations, ...
                                                              Y, 1, size(Y,2), ...
                                                              zeros(mm,1), Pinf, Pstar, ...
                                                              kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                              T,R,Q,H,Z,mm,pp,rr);
        end
        diffuse_periods = length(dlik);
        if isinf(dLIK)
            % Go to univariate diffuse filter if singularity problem.
            singular_diffuse_filter
        end
    end
    if singular_diffuse_filter || (kalman_algo==4)
        % Univariate Diffuse Kalman Filter
        if isequal(H,0)
            H1 = zeros(nobs,1);
            mmm = mm;
        else
            if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                H1 = diag(H);
                mmm = mm; 
            else
                Z = [Z, eye(pp)];
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blckdiag(Pinf,zeros(pp));
                H1 = zeros(nobs,1);
                mmm   = mm+pp;
            end
        end
        % no need to test again for correlation elements
        correlated_errors_have_been_checked = 1;
        
        [dLIK,dlik,a,Pstar] = univariate_kalman_filter_d(DynareDataset.missing.aindex,...
                                                         DynareDataset.missing.number_of_observations,...
                                                         DynareDataset.missing.no_more_missing_observations, ...
                                                         Y, 1, size(Y,2), ...
                                                         zeros(mmm,1), Pinf, Pstar, ...
                                                         kalman_tol, riccati_tol, DynareOptions.presample, ...
                                                         T,R,Q,H1,Z,mmm,pp,rr);
        diffuse_periods = length(dlik);
    end
  case 4% Start from the solution of the Riccati equation.
    if kalman_algo ~= 2
        kalman_algo = 1;
    end
    if isequal(H,0)
        [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,np,length(Z))));
    else
        [err,Pstar] = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,np,length(Z))),H);
    end
    if err
        disp(['DsgeLikelihood:: I am not able to solve the Riccati equation, so I switch to lik_init=1!']);
        DynareOptions.lik_init = 1;
        Pstar = lyapunov_symm(T,R*Q*R',DynareOptions.qz_criterium,DynareOptions.lyapunov_complex_threshold);
    end
    Pinf  = [];
  otherwise
    error('DsgeLikelihood:: Unknown initialization approach for the Kalman filter!')
end

%------------------------------------------------------------------------------
% 4. Likelihood evaluation
%------------------------------------------------------------------------------

if ((kalman_algo==1) || (kalman_algo==3))% Multivariate Kalman Filter
    if no_missing_data_flag
        [LIK,lik] = kalman_filter(Y,diffuse_periods+1,size(Y,2), ...
                            a,Pstar, ...
                            kalman_tol, riccati_tol, ...
                            DynareOptions.presample, ...
                            T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods);
    else
        [LIK,lik] = missing_observations_kalman_filter(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                               a, Pstar, ...
                                               kalman_tol, DynareOptions.riccati_tol, ...
                                               DynareOptions.presample, ...
                                               T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods);
    end
    if isinf(LIK)
        if kalman_algo == 1
            kalman_algo = 2;
        else
            kalman_algo = 4;
        end
        singularity_flag = 1;
    else
        if DynareOptions.lik_init==3
            LIK = LIK + dLIK;
            lik = [dlik; lik];
        end
    end
end

if (kalman_algo==2) || (kalman_algo==4)
    % Univariate Kalman Filter
    % resetting measurement error covariance matrix when necessary                                                           % 
    if ~correlated_errors_have_been_checked
        if isequal(H,0)
            H = zeros(nobs,1);
            mmm = mm;
        else
            if all(all(abs(H-diag(diag(H)))<1e-14))% ie, the covariance matrix is diagonal...
                H = diag(H);
                mmm = mm;
            else
                Z = [Z, eye(pp)];
                T = blkdiag(T,zeros(pp));
                Q = blkdiag(Q,H);
                R = blkdiag(R,eye(pp));
                Pstar = blkdiag(Pstar,H);
                Pinf  = blckdiag(Pinf,zeros(pp));
                H = zeros(nobs,1);
                mmm   = mm+pp;
            end
        end
    end

    [LIK,lik] = univariate_kalman_filter(DynareDataset.missing.aindex,DynareDataset.missing.number_of_observations,DynareDataset.missing.no_more_missing_observations,Y,diffuse_periods+1,size(Y,2), ...
                                       a,Pstar, ...
                                       DynareOptions.kalman_tol, ...
                                       DynareOptions.riccati_tol, ...
                                       DynareOptions.presample, ...
                                       T,Q,R,H,Z,mmm,pp,rr,diffuse_periods);
    if DynareOptions.lik_init==3
        LIK = LIK+dLIK;
        lik = [dlik; lik];
    end
end

if isnan(LIK)
    info = 45;
    exit_flag = 0;
    return
end
if imag(LIK)~=0
    likelihood = penalty;
else
    likelihood = LIK;
end

% ------------------------------------------------------------------------------
% 5. Adds prior if necessary
% ------------------------------------------------------------------------------
lnprior = priordens(xparam1,BayesInfo.pshape,BayesInfo.p6,BayesInfo.p7,BayesInfo.p3,BayesInfo.p4);
fval    = (likelihood-lnprior);

% Update DynareOptions.kalman_algo.
DynareOptions.kalman_algo = kalman_algo;

% Update the penalty.
penalty = fval;

% Add the prior density at the top of the vector for the density of each observation.
lik=lik(start:end,:);
llik=[-lnprior; lik(:)];