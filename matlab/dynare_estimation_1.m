function dynare_estimation_1(var_list_,dname)
% function dynare_estimation_1(var_list_,dname)
% runs the estimation of the model
%
% INPUTS
%   var_list_:  selected endogenous variables vector
%   dname:      alternative directory name
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2013 Dynare Team
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

global M_ options_ oo_ estim_params_ bayestopt_ dataset_

% Set particle filter flag.
if options_.order > 1
    if options_.particle.status && options_.order==2
        skipline()
        disp('Estimation using a non linear filter!')
        skipline()
        if ~options_.nointeractive && ismember(options_.mode_compute,[1,3,4]) % Known gradient-based optimizers
            disp('You are using a gradient-based mode-finder. Particle filtering introduces discontinuities in the') 
            disp('objective function w.r.t the parameters. Thus, should use a non-gradient based optimizer.')
            fprintf('\nPlease choose a mode-finder:\n')
            fprintf('\t 0 - Continue using gradient-based method (it is most likely that you will no get any sensible result).\n')
            fprintf('\t 6 - Monte Carlo based algorithm\n')
            fprintf('\t 7 - Nelder-Mead simplex based optimization routine (Matlab optimization toolbox required)\n')
            fprintf('\t 8 - Nelder-Mead simplex based optimization routine (Dynare''s implementation)\n')
            fprintf('\t 9 - CMA-ES (Covariance Matrix Adaptation Evolution Strategy) algorithm\n')
            choice = [];
            while isempty(choice)
                choice = input('Please enter your choice: ');
                if isnumeric(choice) && isint(choice) && ismember(choice,[0 6 7 8 9])
                    if choice
                        options_.mode_compute = choice;
                    end
                else
                    fprintf('\nThis is an invalid choice (you have to choose between 0, 6, 7, 8 and 9).\n')
                    choice = [];
                end
            end
        end
    elseif options_.particle.status && options_.order>2
        error(['Non linear filter are not implemented with order ' int2str(options_.order) ' approximation of the model!'])
    elseif ~options_.particle.status && options_.order==2
        error('For estimating the model with a second order approximation using a non linear filter, one should have options_.particle.status=1;')
    else
        error(['Cannot estimate a model with an order ' int2str(options_.order) ' approximation!'])
    end
end

if ~options_.dsge_var
    if options_.particle.status
        objective_function = str2func('non_linear_dsge_likelihood');
    else
        objective_function = str2func('dsge_likelihood');
    end
else
    objective_function = str2func('DsgeVarLikelihood');
end

[dataset_,xparam1, hh, M_, options_, oo_, estim_params_,bayestopt_] = dynare_estimation_init(var_list_, dname, [], M_, options_, oo_, estim_params_, bayestopt_);

% Set sigma_e_is_diagonal flag (needed if the shocks block is not declared in the mod file).
M_.sigma_e_is_diagonal = 1;
if estim_params_.ncx || any(nnz(tril(M_.Sigma_e,-1)))
    M_.sigma_e_is_diagonal = 0;
end

% Set the correlation matrix if necessary.
if ~isequal(estim_params_.ncx,nnz(tril(M_.Sigma_e,-1)))
    M_.Correlation_matrix = diag(1./sqrt(diag(M_.Sigma_e)))*M_.Sigma_e*diag(1./sqrt(diag(M_.Sigma_e)));
    % Remove NaNs appearing because of variances calibrated to zero.
    if any(isnan(M_.Correlation_matrix))
        zero_variance_idx = find(~diag(M_.Sigma_e));
        for i=1:length(zero_variance_idx)
            M_.Correlation_matrix(zero_variance_idx(i),:) = 0;
            M_.Correlation_matrix(:,zero_variance_idx(i)) = 0;
        end
    end
end

% Set the correlation matrix of measurement errors if necessary.
if ~isequal(estim_params_.ncn,nnz(tril(M_.H,-1)))
    M_.Correlation_matrix_ME = diag(1./sqrt(diag(M_.H)))*M_.H*diag(1./sqrt(diag(M_.H)));
    % Remove NaNs appearing because of variances calibrated to zero.
    if any(isnan(M_.Correlation_matrix_ME))
        zero_variance_idx = find(~diag(M_.H));
        for i=1:length(zero_variance_idx)
            M_.Correlation_matrix_ME(zero_variance_idx(i),:) = 0;
            M_.Correlation_matrix_ME(:,zero_variance_idx(i)) = 0;
        end
    end
end

data = dataset_.data;
rawdata = dataset_.rawdata;
data_index = dataset_.missing.aindex;
missing_value = dataset_.missing.state;

% Set number of observations
gend = options_.nobs;
% Set the number of observed variables.
n_varobs = size(options_.varobs,1);
% Get the number of parameters to be estimated.
nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.
nx  = nvx+nvn+ncx+ncn+np; % Total number of parameters to be estimated.
%% Set the names of the priors.
pnames = ['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
%% Set parameters bounds
lb = bayestopt_.lb;
ub = bayestopt_.ub;

dr = oo_.dr;

%% load optimal_mh_scale parameter if previous run was with
%% mode_compute=6
mh_scale_fname = [M_.fname '_optimal_mh_scale_parameter.mat'];
if exist(mh_scale_fname)
    if options_.mode_compute == 0
        tmp = load(mh_scale_fname,'Scale');
        bayestopt_.mh_jscale = tmp.Scale;
        clear tmp;
    else
        % remove the file if mode_compute ~= 0
        delete('mh_scale_fname')
    end
end

if ~isempty(estim_params_)
    set_parameters(xparam1);
end

% compute sample moments if needed (bvar-dsge)
if options_.dsge_var
    if dataset_.missing.state
        error('I cannot estimate a DSGE-VAR model with missing observations!')
    end
    if options_.noconstant
        evalin('base',...
               ['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
                'var_sample_moments(options_.first_obs,' ...
                'options_.first_obs+options_.nobs-1,options_.dsge_varlag,-1,' ...
                'options_.datafile, options_.varobs,options_.xls_sheet,options_.xls_range);'])
    else% The steady state is non zero ==> a constant in the VAR is needed!
        evalin('base',['[mYY,mXY,mYX,mXX,Ydata,Xdata] = ' ...
                       'var_sample_moments(options_.first_obs,' ...
                       'options_.first_obs+options_.nobs-1,options_.dsge_varlag,0,' ...
                       'options_.datafile, options_.varobs,options_.xls_sheet,options_.xls_range);'])
    end
end


oo_ = initial_estimation_checks(objective_function,xparam1,dataset_,M_,estim_params_,options_,bayestopt_,oo_);

if isequal(options_.mode_compute,0) && isempty(options_.mode_file) && options_.mh_posterior_mode_estimation==0
    if options_.smoother == 1
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,gend,data,data_index,missing_value);
        oo_.Smoother.SteadyState = ys;
        oo_.Smoother.TrendCoeffs = trend_coeff;
        if options_.filter_covariance
            oo_.Smoother.Variance = P;
        end
        i_endo = bayestopt_.smoother_saved_var_list;
        if options_.nk ~= 0
            oo_.FilteredVariablesKStepAhead = ...
                aK(options_.filter_step_ahead,i_endo,:);
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
        for i=bayestopt_.smoother_saved_var_list'
            i1 = dr.order_var(bayestopt_.smoother_var_list(i));
            eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ...
                  ' = atT(i,:)'';']);
            if options_.nk > 0
                eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ...
                      ' = squeeze(aK(1,i,:));']);
            end
            eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ' = updated_variables(i,:)'';']);
        end
        for i=1:M_.exo_nbr
            eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
        end
    end
    return
end


%% Estimation of the posterior mode or likelihood mode
if ~isequal(options_.mode_compute,0) && ~options_.mh_posterior_mode_estimation
    switch options_.mode_compute
      case 1
        if exist('OCTAVE_VERSION')
            error('Option mode_compute=1 is not available under Octave')
        elseif ~user_has_matlab_license('optimization_toolbox')
            error('Option mode_compute=1 requires the Optimization Toolbox')
        end

        optim_options = optimset('display','iter','LargeScale','off', ...
                                 'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        if options_.analytic_derivation,
            optim_options = optimset(optim_options,'GradObj','on','TolX',1e-7);
        end
            [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
                fmincon(objective_function,xparam1,[],[],[],[],lb,ub,[],optim_options,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
      case 2
        error('ESTIMATION: mode_compute=2 option (Lester Ingber''s Adaptive Simulated Annealing) is no longer available')
      case 3
        if exist('OCTAVE_VERSION') && ~user_has_octave_forge_package('optim')
            error('Option mode_compute=3 requires the optim package')
        elseif ~exist('OCTAVE_VERSION') && ~user_has_matlab_license('optimization_toolbox')
            error('Option mode_compute=3 requires the Optimization Toolbox')
        end

        optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        if options_.analytic_derivation,
            optim_options = optimset(optim_options,'GradObj','on');
        end
        if ~exist('OCTAVE_VERSION')
            [xparam1,fval,exitflag] = fminunc(objective_function,xparam1,optim_options,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        else
            % Under Octave, use a wrapper, since fminunc() does not have a 4th arg
            func = @(x) objective_function(x, dataset_,options_,M_,estim_params_,bayestopt_,oo_);
            [xparam1,fval,exitflag] = fminunc(func,xparam1,optim_options);
        end

      case 4
        H0 = 1e-4*eye(nx);
        crit = 1e-7;
        nit = 1000;
        verbose = 2;
        if options_.analytic_derivation,
            analytic_grad=1;
        else
            analytic_grad=[];
        end

            [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
                csminwel1(objective_function,xparam1,H0,analytic_grad,crit,nit,options_.gradient_method,options_.gradient_epsilon,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
            disp(sprintf('Objective function at mode: %f',fval))
      case 5
        if isfield(options_,'hess')
            flag = options_.hess;
        else
            flag = 1;
        end
        if isfield(options_,'ftol')
            crit = options_.ftol;
        else
            crit = 1.e-5;
        end
        if options_.analytic_derivation,
            analytic_grad=1;
            ana_deriv = options_.analytic_derivation;
            options_.analytic_derivation = -1;
            crit = 1.e-7;
            flag = 0;
        else
            analytic_grad=0;
        end
        if isfield(options_,'nit')
            nit = options_.nit;
        else
            nit=1000;
        end
        [xparam1,hh,gg,fval,invhess] = newrat(objective_function,xparam1,analytic_grad,crit,nit,flag,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        if options_.analytic_derivation,
            options_.analytic_derivation = ana_deriv;
        end
        parameter_names = bayestopt_.name;
        save([M_.fname '_mode.mat'],'xparam1','hh','gg','fval','invhess','parameter_names');
      case 6
        fval = feval(objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        OldMode = fval;
        if ~exist('MeanPar','var')
            MeanPar = xparam1;
        end
        if ~isempty(hh)
            CovJump = inv(hh);
        else% The covariance matrix is initialized with the prior
            % covariance (a diagonal matrix) %%Except for infinite variances ;-)
            varinit = 'prior';
            if strcmpi(varinit,'prior')
                stdev = bayestopt_.p2;
                indx = find(isinf(stdev));
                stdev(indx) = ones(length(indx),1)*sqrt(10);
                vars = stdev.^2;
                CovJump = diag(vars);
            elseif strcmpi(varinit,'eye')
                vars = ones(length(bayestopt_.p2),1)*0.1;
                CovJump = diag(vars);
            else
                disp('gmhmaxlik :: Error!')
                return
            end
        end
        OldPostVar = CovJump;
        Scale = options_.mh_jscale;
        for i=1:options_.gmhmaxlik.iterations
            if i == 1
                if options_.gmhmaxlik.iterations>1
                    flag = '';
                else
                    flag = 'LastCall';
                end
                [xparam1,PostVar,Scale,PostMean] = ...
                    gmhmaxlik(objective_function,xparam1,[lb ub],options_.gmhmaxlik,Scale,flag,MeanPar,CovJump,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                fval = feval(objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                options_.mh_jscale = Scale;
                mouvement = max(max(abs(PostVar-OldPostVar)));
                skipline()
                disp('========================================================== ')
                disp(['   Change in the covariance matrix = ' num2str(mouvement) '.'])
                disp(['   Mode improvement = ' num2str(abs(OldMode-fval))])
                disp(['   New value of jscale = ' num2str(Scale)])
                disp('========================================================== ')
                OldMode = fval;
            else
                OldPostVar = PostVar;
                if i<options_.gmhmaxlik.iterations
                    flag = '';
                else
                    flag = 'LastCall';
                end
                [xparam1,PostVar,Scale,PostMean] = ...
                    gmhmaxlik(objective_function,xparam1,[lb ub],...
                              options_.gmhmaxlik,Scale,flag,PostMean,PostVar,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                fval = feval(objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                options_.mh_jscale = Scale;
                mouvement = max(max(abs(PostVar-OldPostVar)));
                fval = feval(objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                skipline()
                disp('========================================================== ')
                disp(['   Change in the covariance matrix = ' num2str(mouvement) '.'])
                disp(['   Mode improvement = ' num2str(abs(OldMode-fval))])
                disp(['   New value of jscale = ' num2str(Scale)])
                disp('========================================================== ')
                OldMode = fval;
            end
            hh = inv(PostVar);
            save([M_.fname '_mode.mat'],'xparam1','hh');
            save([M_.fname '_optimal_mh_scale_parameter.mat'],'Scale');
            bayestopt_.jscale = ones(length(xparam1),1)*Scale;
        end
        skipline()
        disp(['Optimal value of the scale parameter = ' num2str(Scale)])
        skipline()
        disp(['Final value of the log posterior (or likelihood): ' num2str(fval)])
        skipline()
        parameter_names = bayestopt_.name;
        save([M_.fname '_mode.mat'],'xparam1','hh','parameter_names');
      case 7
        % Matlab's simplex (Optimization toolbox needed).
        if exist('OCTAVE_VERSION') && ~user_has_octave_forge_package('optim')
            error('Option mode_compute=7 requires the optim package')
        elseif ~exist('OCTAVE_VERSION') && ~user_has_matlab_license('optimization_toolbox')
            error('Option mode_compute=7 requires the Optimization Toolbox')
        end

        optim_options = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options_,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        [xparam1,fval,exitflag] = fminsearch(objective_function,xparam1,optim_options,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
      case 8
        % Dynare implementation of the simplex algorithm.
        [xparam1,fval,exitflag] = simplex_optimization_routine(objective_function,xparam1,options_.simplex,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
      case 9
        H0 = 1e-4*ones(nx,1);
        warning('off','CMAES:NonfinitenessRange');
        warning('off','CMAES:InitialSigma');
        [x, fval, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes(func2str(objective_function),xparam1,H0,options_.cmaes,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        xparam1=BESTEVER.x;
        disp(sprintf('\n Objective function at mode: %f',fval))
      case 10 
         options_.cova_compute = 0 ;
         [xparam1,stdh,lb_95,ub_95,med_param] = online_auxiliary_filter(xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_) ;
      case 101
        myoptions=soptions;
        myoptions(2)=1e-6; %accuracy of argument
        myoptions(3)=1e-6; %accuracy of function (see Solvopt p.29)
        myoptions(5)= 1.0;
        [xparam1,fval]=solvopt(xparam1,objective_function,[],myoptions,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
      case 102
        %simulating annealing
        %        LB=zeros(size(xparam1))-20;
        % UB=zeros(size(xparam1))+20;
        LB = xparam1 - 1;
        UB = xparam1 + 1;
        neps=10;
        %  Set input parameters.
        maxy=0;
        epsilon=1.0e-9;
        rt_=.10;
        t=15.0;
        ns=10;
        nt=10;
        maxevl=100000000;
        idisp =1;
        npar=length(xparam1);

        disp(['size of param',num2str(length(xparam1))])
        c=.1*ones(npar,1);
        %*  Set input values of the input/output parameters.*

        vm=1*ones(npar,1);
        disp(['number of parameters= ' num2str(npar) 'max= '  num2str(maxy) 't=  ' num2str(t)]);
        disp(['rt_=  '  num2str(rt_) 'eps=  '  num2str(epsilon) 'ns=  '  num2str(ns)]);
        disp(['nt=  '  num2str(nt) 'neps= '   num2str(neps) 'maxevl=  '  num2str(maxevl)]);
        %      disp(['iprint=   '   num2str(iprint) 'seed=   '   num2str(seed)]);
        disp '  ';
        disp '  ';
        disp(['starting values(x)  ' num2str(xparam1')]);
        disp(['initial step length(vm)  '  num2str(vm')]);
        disp(['lower bound(lb)', 'initial conditions', 'upper bound(ub)' ]);
        disp([LB xparam1 UB]);
        disp(['c vector   ' num2str(c')]);

        [xparam1, fval, nacc, nfcnev, nobds, ier, t, vm] = sa(objective_function,xparam1,maxy,rt_,epsilon,ns,nt ...
                                                              ,neps,maxevl,LB,UB,c,idisp ,t,vm,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
      case 'prior'
        hh = diag(bayestopt_.p2.^2);
      otherwise
        if ischar(options_.mode_compute)
            [xparam1, fval, retcode ] = feval(options_.mode_compute,objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        else
            error(['dynare_estimation:: mode_compute = ' int2str(options_.mode_compute) ' option is unknown!'])
        end
    end
    if ~isequal(options_.mode_compute,6) && ~isequal(options_.mode_compute,'prior')
        if options_.cova_compute == 1
            if options_.analytic_derivation && strcmp(func2str(objective_function),'dsge_likelihood'),
                ana_deriv = options_.analytic_derivation;
                options_.analytic_derivation = 2;
                [junk1, junk2, hh] = feval(objective_function,xparam1, ...
                    dataset_,options_,M_,estim_params_,bayestopt_,oo_);
                options_.analytic_derivation = ana_deriv;
                
            else
                hh = reshape(hessian(objective_function,xparam1, ...
                    options_.gstep,dataset_,options_,M_,estim_params_,bayestopt_,oo_),nx,nx);
            end
        end
    end
    parameter_names = bayestopt_.name;
    if options_.cova_compute
        save([M_.fname '_mode.mat'],'xparam1','hh','parameter_names');
    else
        save([M_.fname '_mode.mat'],'xparam1','parameter_names');
    end
end

if options_.cova_compute == 0
    hh = [];%NaN(length(xparam1),length(xparam1));
end

if ~options_.mh_posterior_mode_estimation && options_.cova_compute
    try
        chol(hh);
    catch
        skipline()
        disp('POSTERIOR KERNEL OPTIMIZATION PROBLEM!')
        disp(' (minus) the hessian matrix at the "mode" is not positive definite!')
        disp('=> posterior variance of the estimated parameters are not positive.')
        disp('You should  try  to change the initial values of the parameters using')
        disp('the estimated_params_init block, or use another optimization routine.')
        params_at_bound=find(xparam1==ub | xparam1==lb);
        if ~isempty(params_at_bound)
            for ii=1:length(params_at_bound)
            params_at_bound_name{ii,1}=get_the_name(ii,0,M_,estim_params_,options_);
            end
            disp_string=[params_at_bound_name{1,:}];
            for ii=2:size(params_at_bound_name,1)
                disp_string=[disp_string,', ',params_at_bound_name{ii,:}];
            end
            fprintf('\nThe following parameters are at the prior bound: %s\n', disp_string)
            fprintf('Some potential solutions are:\n')
            fprintf('   - Check your model for mistakes.\n')
            fprintf('   - Check whether model and data are consistent (correct observation equation).\n')
            fprintf('   - Shut off prior_trunc.\n')
            fprintf('   - Use a different mode_compute like 6 or 9.\n')
            fprintf('   - Check whether the parameters estimated are identified.\n')
            fprintf('   - Increase the informativeness of the prior.\n')
        end
        warning('The results below are most likely wrong!');
    end
end

if options_.mode_check.status == 1 && ~options_.mh_posterior_mode_estimation
    ana_deriv = options_.analytic_derivation;
    options_.analytic_derivation = 0;
    mode_check(objective_function,xparam1,hh,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
    options_.analytic_derivation = ana_deriv;
end

oo_.posterior.optimization.mode = xparam1;
oo_.posterior.optimization.Variance = [];
if ~options_.mh_posterior_mode_estimation
    if options_.cova_compute
        invhess = inv(hh);
        stdh = sqrt(diag(invhess));
        oo_.posterior.optimization.Variance = invhess;
    end
else
    variances = bayestopt_.p2.*bayestopt_.p2;
    idInf = isinf(variances);
    variances(idInf) = 1;
    invhess = options_.mh_posterior_mode_estimation*diag(variances);
    xparam1 = bayestopt_.p5;
    idNaN = isnan(xparam1);
    xparam1(idNaN) = bayestopt_.p1(idNaN);
end

if ~options_.cova_compute
    stdh = NaN(length(xparam1),1);
end

if any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    %% display results table and store parameter estimates and standard errors in results
    oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,'Posterior','posterior');
    %% Laplace approximation to the marginal log density:
    if options_.cova_compute
        estim_params_nbr = size(xparam1,1);
        scale_factor = -sum(log10(diag(invhess)));
        log_det_invhess = -estim_params_nbr*log(scale_factor)+log(det(scale_factor*invhess));
        likelihood = feval(objective_function,xparam1,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        oo_.MarginalDensity.LaplaceApproximation = .5*estim_params_nbr*log(2*pi) + .5*log_det_invhess - likelihood;
        skipline()
        disp(sprintf('Log data density [Laplace approximation] is %f.',oo_.MarginalDensity.LaplaceApproximation))
        skipline()
    end
elseif ~any(bayestopt_.pshape > 0) && ~options_.mh_posterior_mode_estimation
    oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,'Maximum Likelihood','mle');
end

OutputDirectoryName = CheckPath('Output',M_.dname);

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_params.mat'],'pindx');
end

if (any(bayestopt_.pshape  >0 ) && options_.mh_replic) || ...
        (any(bayestopt_.pshape >0 ) && options_.load_mh_file)  %% not ML estimation
    bounds = prior_bounds(bayestopt_,options_);
    bounds(:,1)=max(bounds(:,1),lb);
    bounds(:,2)=min(bounds(:,2),ub);
    bayestopt_.lb = bounds(:,1);
    bayestopt_.ub = bounds(:,2);
    outside_bound_pars=find(xparam1 < bounds(:,1) | xparam1 > bounds(:,2));
    if ~isempty(outside_bound_pars)
        for ii=1:length(outside_bound_pars)
            outside_bound_par_names{ii,1}=get_the_name(ii,0,M_,estim_params_,options_);
        end
        disp_string=[outside_bound_par_names{1,:}];
        for ii=2:size(outside_bound_par_names,1)
            disp_string=[disp_string,', ',outside_bound_par_names{ii,:}];
        end
        error(['Mode value(s) of ', disp_string ,' are outside parameter bounds. Potentially, you should set prior_trunc=0.'])
    end
    % runs MCMC
    if options_.mh_replic
        if options_.load_mh_file && options_.use_mh_covariance_matrix
            invhess = compute_mh_covariance_matrix;
        end
        ana_deriv = options_.analytic_derivation;
        options_.analytic_derivation = 0;
        if options_.cova_compute
            oo_.MC_record=feval(options_.posterior_sampling_method,objective_function,options_.proposal_distribution,xparam1,invhess,bounds,dataset_,options_,M_,estim_params_,bayestopt_,oo_);
        else
            error('I Cannot start the MCMC because the Hessian of the posterior kernel at the mode was not computed.')
        end
        options_.analytic_derivation = ana_deriv;
    end
    if options_.mh_posterior_mode_estimation
        CutSample(M_, options_, estim_params_);
        return
    else
        if ~options_.nodiagnostic && options_.mh_replic > 2000 && options_.mh_nblck > 1
            McMCDiagnostics(options_, estim_params_, M_);
        end
        %% Here i discard first half of the draws:
        CutSample(M_, options_, estim_params_);
        %% Estimation of the marginal density from the Mh draws:
        if options_.mh_replic
            [marginal,oo_] = marginal_density(M_, options_, estim_params_, oo_);
            oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_);
            if ~options_.nograph
                oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_);
            end
            [oo_.posterior.metropolis.mean,oo_.posterior.metropolis.Variance] ...
                = GetPosteriorMeanVariance(M_,options_.mh_drop);
        else
            load([M_.fname '_results'],'oo_');
        end
        metropolis_draw(1);
        if options_.bayesian_irf
            PosteriorIRF('posterior');
        end
        if options_.moments_varendo
            oo_ = compute_moments_varendo('posterior',options_,M_,oo_,var_list_);
        end
        if options_.smoother || ~isempty(options_.filter_step_ahead) || options_.forecast
            prior_posterior_statistics('posterior',dataset_);
        end
        xparam = get_posterior_parameters('mean');
        M_ = set_all_parameters(xparam,estim_params_,M_);
    end
end

if options_.particle.status
    return
end

if (~((any(bayestopt_.pshape > 0) && options_.mh_replic) || (any(bayestopt_.pshape ...
                                                      > 0) && options_.load_mh_file)) ...
    || ~options_.smoother ) && options_.partial_information == 0  % to be fixed
    %% ML estimation, or posterior mode without metropolis-hastings or metropolis without bayesian smooth variable
    [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,T,R,P,PK,decomp] = DsgeSmoother(xparam1,dataset_.info.ntobs,dataset_.data,dataset_.missing.aindex,dataset_.missing.state);
    oo_.Smoother.SteadyState = ys;
    oo_.Smoother.TrendCoeffs = trend_coeff;
    oo_.Smoother.Variance = P;
    i_endo = bayestopt_.smoother_saved_var_list;
    if options_.nk ~= 0
        oo_.FilteredVariablesKStepAhead = aK(options_.filter_step_ahead, ...
                                             i_endo,:);
        if isfield(options_,'kalman_algo')
            if ~isempty(PK)
                oo_.FilteredVariablesKStepAheadVariances = ...
                    PK(options_.filter_step_ahead,i_endo,i_endo,:);
            end
            if ~isempty(decomp)
                oo_.FilteredVariablesShockDecomposition = ...
                    decomp(options_.filter_step_ahead,i_endo,:,:);
            end
        end
    end
    for i=bayestopt_.smoother_saved_var_list'
        i1 = dr.order_var(bayestopt_.smoother_var_list(i));
        eval(['oo_.SmoothedVariables.' deblank(M_.endo_names(i1,:)) ' = ' ...
                            'atT(i,:)'';']);
        if options_.nk > 0
            eval(['oo_.FilteredVariables.' deblank(M_.endo_names(i1,:)) ...
                  ' = squeeze(aK(1,i,:));']);
        end
        eval(['oo_.UpdatedVariables.' deblank(M_.endo_names(i1,:)) ...
              ' = updated_variables(i,:)'';']);
    end
    for i=1:M_.exo_nbr
        eval(['oo_.SmoothedShocks.' deblank(M_.exo_names(i,:)) ' = innov(i,:)'';']);
    end
    if ~options_.nograph,
        [nbplt,nr,nc,lr,lc,nstar] = pltorg(M_.exo_nbr);
        if options_.TeX
            fidTeX = fopen([M_.fname '_SmoothedShocks.TeX'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        for plt = 1:nbplt,
            fh = dyn_figure(options_,'Name','Smoothed shocks');
            NAMES = [];
            if options_.TeX, TeXNAMES = []; end
            nstar0=min(nstar,M_.exo_nbr-(plt-1)*nstar);
            if gend==1
                marker_string{1,1}='-ro';
                marker_string{2,1}='-ko';
            else
                marker_string{1,1}='-r';
                marker_string{2,1}='-k';
            end
            for i=1:nstar0,
                k = (plt-1)*nstar+i;
                subplot(nr,nc,i);
                plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                hold on
                plot(1:gend,innov(k,:),marker_string{2,1},'linewidth',1)
                hold off
                name = deblank(M_.exo_names(k,:));
                if isempty(NAMES)
                    NAMES = name;
                else
                    NAMES = char(NAMES,name);
                end
                if ~isempty(options_.XTick)
                    set(gca,'XTick',options_.XTick)
                    set(gca,'XTickLabel',options_.XTickLabel)
                end
                if gend>1
                    xlim([1 gend])
                end
                if options_.TeX
                    texname = M_.exo_names_tex(k,:);
                    if isempty(TeXNAMES)
                        TeXNAMES = ['$ ' deblank(texname) ' $'];
                    else
                        TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                    end
                end
                title(name,'Interpreter','none')
            end
            dyn_saveas(fh,[M_.fname '_SmoothedShocks' int2str(plt)],options_);
            if options_.TeX
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                for jj = 1:nstar0
                    fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                end
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',M_.fname,int2str(plt));
                fprintf(fidTeX,'\\caption{Smoothed shocks.}');
                fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,'\n');
            end
        end
        if options_.TeX
            fprintf(fidTeX,'\n');
            fprintf(fidTeX,'%% End of TeX file.\n');
            fclose(fidTeX);
        end
    end
    %%
    %%  Smooth observational errors...
    %%
    if options_.prefilter == 1
        yf = atT(bayestopt_.mf,:)+repmat(dataset_.descriptive.mean',1,gend);
    elseif options_.loglinear == 1
        yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
             trend_coeff*[1:gend];
    else
        yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
             trend_coeff*[1:gend];
    end
    if nvn
        number_of_plots_to_draw = 0;
        index = [];
        for i=1:n_varobs
            if max(abs(measurement_error)) > 0.000000001
                number_of_plots_to_draw = number_of_plots_to_draw + 1;
                index = cat(1,index,i);
            end
            eval(['oo_.SmoothedMeasurementErrors.' deblank(options_.varobs(i,:)) ...
                  ' = measurement_error(i,:)'';']);
        end
        if ~options_.nograph
            [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
            if options_.TeX
                fidTeX = fopen([M_.fname '_SmoothedObservationErrors.TeX'],'w');
                fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
                fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
                fprintf(fidTeX,' \n');
            end
            for plt = 1:nbplt
                fh = dyn_figure(options_,'Name','Smoothed observation errors');
                NAMES = [];
                if options_.TeX, TeXNAMES = []; end
                nstar0=min(nstar,number_of_plots_to_draw-(nbplt-1)*nstar);
                if gend==1
                    marker_string{1,1}='-ro';
                    marker_string{2,1}='-ko';
                else
                    marker_string{1,1}='-r';
                    marker_string{2,1}='-k';
                end
                for i=1:nstar0
                    k = (plt-1)*nstar+i;
                    subplot(nr,nc,i);
                    plot([1 gend],[0 0],marker_string{1,1},'linewidth',.5)
                    hold on
                    plot(1:gend,measurement_error(index(k),:),marker_string{2,1},'linewidth',1)
                    hold off
                    name = deblank(options_.varobs(index(k),:));
                    if gend>1
                        xlim([1 gend])
                    end
                    if isempty(NAMES)
                        NAMES = name;
                    else
                        NAMES = char(NAMES,name);
                    end
                    if ~isempty(options_.XTick)
                        set(gca,'XTick',options_.XTick)
                        set(gca,'XTickLabel',options_.XTickLabel)
                    end
                    if options_.TeX
                        idx = strmatch(options_.varobs(index(k),:),M_.endo_names,'exact');
                        texname = M_.endo_names_tex(idx,:);
                        if isempty(TeXNAMES)
                            TeXNAMES = ['$ ' deblank(texname) ' $'];
                        else
                            TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                        end
                    end
                    title(name,'Interpreter','none')
                end
                dyn_saveas(fh,[M_.fname '_SmoothedObservationErrors' int2str(plt)],options_);
                if options_.TeX
                    fprintf(fidTeX,'\\begin{figure}[H]\n');
                    for jj = 1:nstar0
                        fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
                    end
                    fprintf(fidTeX,'\\centering \n');
                    fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',M_.fname,int2str(plt));
                    fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
                    fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
                    fprintf(fidTeX,'\\end{figure}\n');
                    fprintf(fidTeX,'\n');
                end
            end
            if options_.TeX
                fprintf(fidTeX,'\n');
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end
        end
    end
    %%
    %%  Historical and smoothed variabes
    %%
    if ~options_.nograph
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
    if options_.TeX
        fidTeX = fopen([M_.fname '_HistoricalAndSmoothedVariables.TeX'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
    end
    for plt = 1:nbplt,
        fh = dyn_figure(options_,'Name','Historical and smoothed variables');
        NAMES = [];
        if options_.TeX, TeXNAMES = []; end
        nstar0=min(nstar,n_varobs-(plt-1)*nstar);
        if gend==1
           marker_string{1,1}='-ro';
           marker_string{2,1}='--ko';
        else
           marker_string{1,1}='-r';
           marker_string{2,1}='--k';
        end
        for i=1:nstar0,
            k = (plt-1)*nstar+i;
            subplot(nr,nc,i);
            plot(1:gend,yf(k,:),marker_string{1,1},'linewidth',1)
            hold on
            plot(1:gend,rawdata(:,k),marker_string{2,1},'linewidth',1)
            hold off
            name = deblank(options_.varobs(k,:));
            if isempty(NAMES)
                NAMES = name;
            else
                NAMES = char(NAMES,name);
            end
            if ~isempty(options_.XTick)
                set(gca,'XTick',options_.XTick)
                set(gca,'XTickLabel',options_.XTickLabel)
            end
            if gend>1
                xlim([1 gend])
            end
            if options_.TeX
                idx = strmatch(options_.varobs(k,:),M_.endo_names,'exact');
                texname = M_.endo_names_tex(idx,:);
                if isempty(TeXNAMES)
                    TeXNAMES = ['$ ' deblank(texname) ' $'];
                else
                    TeXNAMES = char(TeXNAMES,['$ ' deblank(texname) ' $']);
                end
            end
            title(name,'Interpreter','none')
        end
        dyn_saveas(fh,[M_.fname '_HistoricalAndSmoothedVariables' int2str(plt)],options_);
        if options_.TeX
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for jj = 1:nstar0,
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
            end
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',M_.fname,int2str(plt));
            fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
            fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,'\n');
        end
    end
    if options_.TeX
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    end
end

if options_.forecast > 0 && options_.mh_replic == 0 && ~options_.load_mh_file
    dyn_forecast(var_list_,'smoother');
end

if np > 0
    pindx = estim_params_.param_vals(:,1);
    save([M_.fname '_pindx.mat'] ,'pindx');
end

