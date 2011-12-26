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

% Copyright (C) 2009, 2010, 2011 Dynare Team
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
    
debug = 0;
options_.verbosity = options_.ep.verbosity;
verbosity = options_.ep.verbosity+debug;

% Test if bytecode and block options are used (these options are mandatory)
if ~( options_.bytecode && options_.block )
    error('extended_path:: Options bytecode and block are mandatory!')
end

% Set default initial conditions.
if isempty(initial_conditions)
    initial_conditions = oo_.steady_state;
end

% Set maximum number of iterations for the deterministic solver.
options_.maxit_ = options_.ep.maxit;

% Set the number of periods for the perfect foresight model
options_.periods = options_.ep.periods;

% Set the algorithm for the perfect foresight solver
options_.stack_solve_algo = options_.ep.stack_solve_algo;

% Compute the first order reduced form if needed.
%
% REMARK. It is assumed that the user did run the same mod file with stoch_simul(order=1) and save
% all the globals in a mat file called linear_reduced_form.mat;
if options_.ep.init
    lrf = load('linear_reduced_form','oo_');
    oo_.dr = lrf.oo_.dr; clear('lrf');
    if options_.ep.init==2
        lambda = .8;
    end
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
options_.minimal_solving_period = options_.ep.periods;

% Get indices of variables with non zero steady state
idx = find(abs(oo_.steady_state)>0);

% Initialize the exogenous variables.
make_ex_;

% Initialize the endogenous variables.
make_y_;

% Initialize the output array.
time_series = zeros(M_.endo_nbr,sample_size);

% Set the covariance matrix of the structural innovations.
variances = diag(M_.Sigma_e); 
positive_var_indx = find(variances>0); 
covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx); 
number_of_structural_innovations = length(covariance_matrix); 
covariance_matrix_upper_cholesky = chol(covariance_matrix); 

% Set seed.
if options_.ep.set_dynare_seed_to_default
    set_dynare_seed('default');
end

% Simulate shocks.
switch options_.ep.innovation_distribution
  case 'gaussian'
      oo_.ep.shocks = randn(sample_size,number_of_structural_innovations)*covariance_matrix_upper_cholesky; 
  otherwise
    error(['extended_path:: ' options_.ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
end
    
% Set future shocks (Stochastic Extended Path approach)
if options_.ep.stochastic.status
    switch options_.ep.stochastic.method
      case 'tensor'
        [r,w] = gauss_hermite_weights_and_nodes(options_.ep.stochastic.nodes);
        if options_.ep.stochastic.order*M_.exo_nbr>1
            for i=1:options_.ep.stochastic.order*M_.exo_nbr
                rr(i) = {r};
                ww(i) = {w};
            end
            rrr = cartesian_product_of_sets(rr{:});
            www = cartesian_product_of_sets(ww{:});
        else
            rrr = r;
            www = w;
        end
        www = prod(www,2);
        number_of_nodes = length(www);
        relative_weights = www/max(www);
        switch options_.ep.stochastic.pruned.status
          case 1
            jdx = find(relative_weights>options_.ep.stochastic.pruned.relative);
            www = www(jdx);
            www = www/sum(www);
            rrr = rrr(jdx,:);
          case 2
            jdx = find(weights>options_.ep.stochastic.pruned.level);
            www = www(jdx);
            www = www/sum(www);
            rrr = rrr(jdx,:);
          otherwise
            % Nothing to be done!
        end
        nnn = length(www);
      otherwise
        error('extended_path:: Unknown stochastic_method option!')
    end
else
    rrr = zeros(1,number_of_structural_innovations);
    www = 1;
    nnn = 1;
end

% Initializes some variables.
t  = 0;

% Set waitbar (graphic or text  mode)
graphic_waitbar_flag = ~( options_.console_mode || exist('OCTAVE_VERSION') );

if graphic_waitbar_flag
    hh = waitbar(0,['Please wait. Extended Path simulations...']);
    set(hh,'Name','EP simulations');
else
    for i=1:2, disp(' '), end
    if ~exist('OCTAVE_VERSION')
       back = [];
    end
end


% Main loop.
while (t<sample_size)
    if ~mod(t,10)
        if graphic_waitbar_flag
            waitbar(t/sample_size);
        else
            if exist('OCTAVE_VERSION')
                printf('Please wait. Extended Path simulations... %3.f%%\r done', 100*t/sample_size);
            else
                str = sprintf('Please wait. Extended Path simulations... %3.f%% done.', 100*t/sample_size);
                fprintf([back '%s'],str);
                back=repmat('\b',1,length(str));
            end
        end
    end
    % Set period index.
    t = t+1;
    shocks = oo_.ep.shocks(t,:);
    % Put it in oo_.exo_simul (second line).
    oo_.exo_simul(2,positive_var_indx) = shocks;
    for s = 1:nnn
        for u=1:options_.ep.stochastic.order
            oo_.exo_simul(2+u,positive_var_indx) = rrr(s,(((u-1)*M_.exo_nbr)+1):(u*M_.exo_nbr))*covariance_matrix_upper_cholesky;
        end
        oo_.exo_simul(2:5,:)
        if options_.ep.init% Compute first order solution... 
            initial_path = simult_(initial_conditions,oo_.dr,oo_.exo_simul(2:end,:),1);
            if options_.ep.init==1
                oo_.endo_simul(:,1:end-1) = initial_path(:,1:end-1);% Last column is the steady state.
            elseif options_.ep.init==2
                oo_.endo_simul(:,1:end-1) = initial_path(:,1:end-1)*lambda+oo_.endo_simul(:,1:end-1)*(1-lambda);
            end
        end
        % Solve a perfect foresight model (using bytecoded version).
        increase_periods = 0;
        endo_simul = oo_.endo_simul;
        while 1
            if ~increase_periods
                t0 = tic;
                [flag,tmp] = bytecode('dynamic'); 
                ctime = toc(t0);
                info.convergence = ~flag;
                info.time = ctime;
            end
            if verbosity
                if info.convergence
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
            % Test if periods is big enough.
            if ~increase_periods &&  max(max(abs(tmp(idx,end-options_.ep.lp:end)./tmp(idx,end-options_.ep.lp-1:end-1)-1)))<options_.dynatol.x
                break
            else
                options_.periods = options_.periods + options_.ep.step;
                options_.minimal_solving_period = options_.periods;
                increase_periods = increase_periods + 1;
                if verbosity
                    if t<10
                        disp(['Time:    ' int2str(t)  '. I increase the number of periods to ' int2str(options_.periods) '.'])
                    elseif t<100
                        disp(['Time:   ' int2str(t) '. I increase the number of periods to ' int2str(options_.periods) '.'])
                    elseif t<1000
                        disp(['Time:  ' int2str(t)  '. I increase the number of periods to ' int2str(options_.periods) '.'])
                    else
                        disp(['Time: ' int2str(t)  '. I increase the number of periods to ' int2str(options_.periods) '.'])
                    end
                end
                if info.convergence
                    oo_.endo_simul = [ tmp , repmat(oo_.steady_state,1,options_.ep.step) ];
                    oo_.exo_simul  = [ oo_.exo_simul ; zeros(options_.ep.step,size(shocks,2)) ];
                    tmp_old = tmp;
                else
                    oo_.endo_simul = [ oo_.endo_simul , repmat(oo_.steady_state,1,options_.ep.step) ];
                    oo_.exo_simul  = [ oo_.exo_simul ; zeros(options_.ep.step,size(shocks,2)) ];
                end
                t0 = tic;
                [flag,tmp] = bytecode('dynamic');
                ctime = toc(t0);
                info.time = info.time+ctime;
                if info.convergence
                    maxdiff = max(max(abs(tmp(:,2:options_.ep.fp)-tmp_old(:,2:options_.ep.fp))));
                    if maxdiff<options_.dynatol.x
                        options_.periods = options_.ep.periods;
                        options_.minimal_solving_period = options_.periods;
                        oo_.exo_simul = oo_.exo_simul(1:(options_.periods+2),:);
                        break
                    end
                else
                    info.convergence = ~flag;
                    if info.convergence
                        continue
                    else
                        if increase_periods==10;
                            if verbosity
                                if t<10
                                    disp(['Time:    ' int2str(t)  '. Even with ' int2str(options_.periods) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                                elseif t<100
                                    disp(['Time:   ' int2str(t)  '. Even with ' int2str(options_.periods) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                                elseif t<1000
                                    disp(['Time:  ' int2str(t)  '. Even with ' int2str(options_.periods) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                                else
                                    disp(['Time: ' int2str(t)  '. Even with ' int2str(options_.periods) ', I am not able to solve the perfect foresight model. Use homotopy instead...'])
                                end
                            end
                            break
                        end
                    end
                end
            end
        end
        if ~info.convergence% If the previous step was unsuccesfull, use an homotopic approach
            [INFO,tmp] = homotopic_steps(.5,.01,t);
            % Cumulate time.
            info.time = ctime+INFO.time;
            if (~isstruct(INFO) && isnan(INFO)) || ~INFO.convergence
                disp('Homotopy:: No convergence of the perfect foresight model solver!')
                error('I am not able to simulate this model!');
            else
                info.convergence = 1;
                oo_.endo_simul = tmp;
                if verbosity && info.convergence
                    disp('Homotopy:: Convergence of the perfect foresight model solver!')
                end
            end
        else
            oo_.endo_simul = tmp;
        end
        % Save results of the perfect foresight model solver.
        time_series(:,t) = time_series(:,t)+ www(s)*oo_.endo_simul(:,2);
        %save('simulated_paths.mat','time_series');
        % Set initial condition for the nex round.
        %initial_conditions = oo_.endo_simul(:,2);
    end
    %oo_.endo_simul = oo_.endo_simul(:,1:options_.periods+M_.maximum_endo_lag+M_.maximum_endo_lead);
    oo_.endo_simul(:,1:end-1) = oo_.endo_simul(:,2:end);
    oo_.endo_simul(:,1) = time_series(:,t);
    oo_.endo_simul(:,end) = oo_.steady_state;
end

if graphic_waitbar_flag
    close(hh);
else
    if ~exist('OCTAVE_VERSION')
        fprintf(back);
    end
end

oo_.endo_simul = oo_.steady_state;