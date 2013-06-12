function time_series = extended_path(initial_conditions,sample_size)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
%  o initial_conditions     [double]    m*nlags array, where m is the number of endogenous variables in the model and
%                                       nlags is the maximum number of lags.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%
% OUTPUTS
%  o time_series            [double]    m*sample_size array, the simulations.
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2009-2013 Dynare Team
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
global M_ options_ oo_

options_.verbosity = options_.ep.verbosity;
verbosity = options_.ep.verbosity+options_.ep.debug;

% Prepare a structure needed by the matlab implementation of the perfect foresight model solver
pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_,'Tensor-Gaussian-Quadrature');

exo_nbr = M_.exo_nbr;
periods = options_.periods;
ep = options_.ep;
steady_state = oo_.steady_state;
dynatol = options_.dynatol;

% Set default initial conditions.
if isempty(initial_conditions)
    initial_conditions = oo_.steady_state;
end

% Set maximum number of iterations for the deterministic solver.
options_.maxit_ = options_.ep.maxit;

% Set the number of periods for the perfect foresight model
periods = options_.ep.periods;
pfm.periods = options_.ep.periods;
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);

% keep a copy of pfm.i_upd
i_upd = pfm.i_upd;

% Set the algorithm for the perfect foresight solver
options_.stack_solve_algo = options_.ep.stack_solve_algo;

% Set check_stability flag
do_not_check_stability_flag = ~options_.ep.check_stability;

% Compute the first order reduced form if needed.
%
% REMARK. It is assumed that the user did run the same mod file with stoch_simul(order=1) and save
% all the globals in a mat file called linear_reduced_form.mat;

dr = struct();
if options_.ep.init
    options_.order = 1;
    [dr,Info,M_,options_,oo_] = resol(1,M_,options_,oo_);
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
options_.minimal_solving_period = 100;%options_.ep.periods;

% Initialize the exogenous variables.
make_ex_;

% Initialize the endogenous variables.
make_y_;

% Initialize the output array.
time_series = zeros(M_.endo_nbr,sample_size);

% Set the covariance matrix of the structural innovations.
variances = diag(M_.Sigma_e);
positive_var_indx = find(variances>0);
effective_number_of_shocks = length(positive_var_indx);
stdd = sqrt(variances(positive_var_indx));
covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx);
covariance_matrix_upper_cholesky = chol(covariance_matrix);

% (re)Set exo_nbr
%exo_nbr = effective_number_of_shocks;

% Set seed.
if options_.ep.set_dynare_seed_to_default
    set_dynare_seed('default');
end

% Set bytecode flag
bytecode_flag = options_.ep.use_bytecode;

% Simulate shocks.
switch options_.ep.innovation_distribution
  case 'gaussian'
      oo_.ep.shocks = randn(sample_size,effective_number_of_shocks)*covariance_matrix_upper_cholesky;
  otherwise
    error(['extended_path:: ' options_.ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
end

% Initializes some variables.
t  = 0;

% Set waitbar (graphic or text  mode)
hh = dyn_waitbar(0,'Please wait. Extended Path simulations...');
set(hh,'Name','EP simulations.');

% hybrid correction
pfm.hybrid_order = options_.ep.stochastic.hybrid_order;
if pfm.hybrid_order
    oo_.dr = set_state_space(oo_.dr,M_,options_);
    options = options_;
    options.order = pfm.hybrid_order;
    pfm.dr = resol(0,M_,options,oo_);
else
    pfm.dr = [];
end

% Main loop.
while (t<sample_size)
    if ~mod(t,10)
        dyn_waitbar(t/sample_size,hh,'Please wait. Extended Path simulations...');
    end
    % Set period index.
    t = t+1;
    shocks = oo_.ep.shocks(t,:);
    % Put it in oo_.exo_simul (second line).
    oo_.exo_simul(2,positive_var_indx) = shocks;
    periods1 = periods;
    exo_simul_1 = zeros(periods1+2,exo_nbr);
    exo_simul_1(2,:) = oo_.exo_simul(2,:);
    pfm1 = pfm;
    info_convergence = 0;
    if ep.init% Compute first order solution (Perturbation)...
        ex = zeros(size(endo_simul_1,2),size(exo_simul_1,2));
        ex(1:size(exo_simul_1,1),:) = exo_simul_1;
        exo_simul_1 = ex;
        initial_path = simult_(initial_conditions,dr,exo_simul_1(2:end,:),1);
        endo_simul_1(:,1:end-1) = initial_path(:,1:end-1)*ep.init+endo_simul_1(:,1:end-1)*(1-ep.init);
    else
        if t==1
            endo_simul_1 = repmat(steady_state,1,periods1+2);
        end
    end
    % Solve a perfect foresight model.
    increase_periods = 0;
    % Keep a copy of endo_simul_1
    endo_simul = endo_simul_1;
    while 1
        if ~increase_periods
            if bytecode_flag && ~options_.ep.stochastic.order
                [flag,tmp] = bytecode('dynamic',endo_simul_1,exo_simul_1);
            else
                flag = 1;
            end
            if flag
                if options_.ep.stochastic.order == 0
                    [flag,tmp,err] = solve_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm1);
                else
                    switch(options_.ep.stochastic.algo)
                        case 0
                        [flag,tmp] = ...
                            solve_stochastic_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm1,options_.ep.stochastic.quadrature.nodes,options_.ep.stochastic.order);
                        case 1
                          [flag,tmp] = ...
                              solve_stochastic_perfect_foresight_model_1(endo_simul_1,exo_simul_1,pfm1,options_.ep.stochastic.quadrature.nodes,options_.ep.stochastic.order);
                    end
                end
            end
            info_convergence = ~flag;
        end
        if verbosity
            if info_convergence
                if t<10
                    disp(['Time:    ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
                elseif t<100
                    disp(['Time:   ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
                elseif t<1000
                    disp(['Time:  ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
                else
                    disp(['Time: ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
                end
            else
                if t<10
                    disp(['Time:    ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
                elseif t<100
                    disp(['Time:   ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
                elseif t<1000
                    disp(['Time:  ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
                else
                    disp(['Time: ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
                end
            end
        end
        if do_not_check_stability_flag
            % Exit from the while loop.
            endo_simul_1 = tmp;
            break
        else
            % Test if periods is big enough.
            % Increase the number of periods.
            periods1 = periods1 + ep.step;
            pfm1.periods = periods1;
            pfm1.i_upd = pfm1.ny+(1:pfm1.periods*pfm1.ny);
            % Increment the counter.
            increase_periods = increase_periods + 1;
            if verbosity
                if t<10
                    disp(['Time:    ' int2str(t)  '. I increase the number of periods to ' int2str(periods1) '.'])
                elseif t<100
                    disp(['Time:   ' int2str(t) '. I increase the number of periods to ' int2str(periods1) '.'])
                elseif t<1000
                    disp(['Time:  ' int2str(t)  '. I increase the number of periods to ' int2str(periods1) '.'])
                else
                    disp(['Time: ' int2str(t)  '. I increase the number of periods to ' int2str(periods1) '.'])
                end
            end
            if info_convergence
                % If the previous call to the perfect foresight model solver exited
                % announcing that the routine converged, adapt the size of endo_simul_1
                % and exo_simul_1.
                endo_simul_1 = [ tmp , repmat(steady_state,1,ep.step) ];
                exo_simul_1  = [ exo_simul_1 ; zeros(ep.step,exo_nbr)];
                tmp_old = tmp;
            else
                % If the previous call to the perfect foresight model solver exited
                % announcing that the routine did not converge, then tmp=1... Maybe
                % should change that, because in some circonstances it may usefull
                % to know where the routine did stop, even if convergence was not
                % achieved.
                endo_simul_1 = [ endo_simul_1 , repmat(steady_state,1,ep.step) ];
                exo_simul_1  = [ exo_simul_1 ; zeros(ep.step,exo_nbr)];
            end
            % Solve the perfect foresight model with an increased number of periods.
            if bytecode_flag && ~options_.ep.stochastic.order
                [flag,tmp] = bytecode('dynamic',endo_simul_1,exo_simul_1);
            else
                flag = 1;
            end
            if flag
                if options_.ep.stochastic.order == 0
                    [flag,tmp,err] = solve_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm1);
                else
                    [flag,tmp] = solve_stochastic_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm1,options_.ep.stochastic.nodes,options_.ep.stochastic.order);
                end
            end
            info_convergence = ~flag;
            if info_convergence
                % If the solver achieved convergence, check that simulated paths did not
                % change during the first periods.
                % Compute the maximum deviation between old path and new path over the
                % first periods
                delta = max(max(abs(tmp(:,2)-tmp_old(:,2))));
                if delta < dynatol.x
                    % If the maximum deviation is close enough to zero, reset the number
                    % of periods to ep.periods
                    periods1 = ep.periods;
                    pfm1.periods = periods1;
                    pfm1.i_upd = pfm1.ny+(1:pfm1.periods*pfm1.ny);
                    % Cut exo_simul_1 and endo_simul_1 consistently with the resetted
                    % number of periods and exit from the while loop.
                    exo_simul_1 = exo_simul_1(1:(periods1+2),:);
                    endo_simul_1 = endo_simul_1(:,1:(periods1+2));
                    break
                end
            else
                % The solver did not converge... Try to solve the model again with a bigger
                % number of periods, except if the number of periods has been increased more
                % than 10 times.
                if increase_periods==10;
                    if verbosity
                        if t<10
                            disp(['Time:    ' int2str(t)  '. Even with ' int2str(periods1) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                        elseif t<100
                            disp(['Time:   ' int2str(t)  '. Even with ' int2str(periods1) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                        elseif t<1000
                            disp(['Time:  ' int2str(t)  '. Even with ' int2str(periods1) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                        else
                            disp(['Time: ' int2str(t)  '. Even with ' int2str(periods1) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                        end
                    end
                    % Exit from the while loop.
                    break
                end
            end% if info_convergence
        end
    end% while
    if ~info_convergence% If exited from the while loop without achieving convergence, use an homotopic approach
        if ~do_not_check_stability_flag
            periods1 = ep.periods;
            pfm1.periods = periods1;
            pfm1.i_upd = i_upd;
            exo_simul_1 = exo_simul_1(1:(periods1+2),:);
            endo_simul_1 = endo_simul_1(:,1:(periods1+2));
        end
        [INFO,tmp] = homotopic_steps(endo_simul,exo_simul_1,.5,.01,pfm1);
        if isstruct(INFO)
            info_convergence = INFO.convergence;
        else
            info_convergence = 0;
        end
        if ~info_convergence
            [INFO,tmp] = homotopic_steps(endo_simul,exo_simul_1,0,.01,pfm1);
            if isstruct(INFO)
                info_convergence = INFO.convergence;
            else
                info_convergence = 0;
            end
            if ~info_convergence
                disp('Homotopy:: No convergence of the perfect foresight model solver!')
                error('I am not able to simulate this model!');
            else
                endo_simul_1 = tmp;
                if verbosity && info_convergence
                    disp('Homotopy:: Convergence of the perfect foresight model solver!')
                end
            end
        else
            info_convergence = 1;
            endo_simul_1 = tmp;
            if verbosity && info_convergence
                disp('Homotopy:: Convergence of the perfect foresight model solver!')
            end
        end
    end
    % Save results of the perfect foresight model solver.
    time_series(:,t) = endo_simul_1(:,2);
    endo_simul_1(:,1:end-1) = endo_simul_1(:,2:end);
    endo_simul_1(:,1) = time_series(:,t);
    endo_simul_1(:,end) = oo_.steady_state;
end% (while) loop over t

dyn_waitbar_close(hh);

oo_.endo_simul = oo_.steady_state;
