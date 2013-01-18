function det_cond_forecast(constrained_paths, constrained_vars, options_cond_fcst, constrained_perfect_foresight)
% Computes conditional forecasts for a deterministic model.
%
% INPUTS
%  o constrained_paths    [double]      m*p array, where m is the number of constrained endogenous variables and p is the number of constrained periods.
%  o constrained_vars     [char]        m*x array holding the names of the controlled endogenous variables.
%  o options_cond_fcst    [structure]   containing the options. The fields are:
%                                                             + replic              [integer]   scalar, number of monte carlo simulations.
%                                                             + parameter_set       [char]      values of the estimated parameters:
%                                                                                               "posterior_mode",
%                                                                                               "posterior_mean",
%                                                                                               "posterior_median",
%                                                                                               "prior_mode" or
%                                                                                               "prior mean".
%                                                                                   [double]     np*1 array, values of the estimated parameters.
%                                                             + controlled_varexo   [char]       m*x array, list of controlled exogenous variables.
%                                                             + conf_sig            [double]     scalar in [0,1], probability mass covered by the confidence bands.
%  o constrained_perfect_foresight [double] m*1 array indicating if the endogenous variables path is perfectly foresight (1) or is a surprise (0)
%
%
% OUTPUTS
%  None.
%
% SPECIAL REQUIREMENTS
%  This routine has to be called after an estimation statement or an estimated_params block.
%
% REMARKS
%  [1] Results are stored in a structure which is saved in a mat file called conditional_forecasts.mat.
%  [2] Use the function plot_icforecast to plot the results.

% Copyright (C) 2013 Dynare Team
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

global options_ oo_ M_

if ~isfield(options_cond_fcst,'periods') || isempty(options_cond_fcst.periods)
    options_cond_fcst.periods = 100;
end

maximum_lag = M_.maximum_lag;
maximum_lead = M_.maximum_lead;
ys = oo_.steady_state;
ny = size(ys,1);
xs = [oo_.exo_steady_state ; oo_.exo_det_steady_state];
nx = size(xs,1);

constrained_periods = size(constrained_paths,2);
n_endo_constrained = size(constrained_vars,1);
if isfield(options_cond_fcst,'controlled_varexo')
    n_control_exo = size(options_cond_fcst.controlled_varexo, 1);
    if n_control_exo ~= n_endo_constrained
        error(['det_cond_forecast:: the number of exogenous controlled variables (' int2str(n_control_exo) ') has to be equal to the number of constrained endogenous variabes (' int2str(n_endo_constrained) ')'])
    end;
else
    error('det_cond_forecast:: to run a deterministic conditional forecast you have to specified the exogenous variables controlled using the option controlled_varex in forecast command');
end;

exo_names = M_.exo_names;
controlled_varexo = zeros(1,n_control_exo);
for i = 1:nx
    for j=1:n_control_exo
        if strcmp(deblank(exo_names(i,:)), deblank(options_cond_fcst.controlled_varexo(j,:)))
            controlled_varexo(j) = i;
        end
    end
end

%todo check if zero => error message

save_options_initval_file = options_.initval_file;
options_.initval_file = '__';

[pos_constrained_pf, junk] = find(constrained_perfect_foresight);
indx_endo_solve_pf = constrained_vars(pos_constrained_pf);
if isempty(indx_endo_solve_pf)
    pf = 0;
else
    pf = length(indx_endo_solve_pf);
end;
indx_endo_solve_surprise = setdiff(constrained_vars, indx_endo_solve_pf);

if isempty(indx_endo_solve_surprise)
    surprise = 0;
else
    surprise = length(indx_endo_solve_surprise);
end;

eps = options_.solve_tolf;
maxit = options_.solve_maxit;

initial_conditions = oo_.steady_state;
past_val = 0;
save_options_periods = options_.periods;
options_.periods = options_cond_fcst.periods;
save_options_dynatol_f = options_.dynatol.f;
options_.dynatol.f = 1e-8;
eps1  = 1e-4;
exo = zeros(maximum_lag + options_cond_fcst.periods, nx);
endo = zeros(maximum_lag + options_cond_fcst.periods, ny);
endo(1,:) = oo_.steady_state';
% if all the endogenous paths are perfectly anticipated we do not need to
% implement extended path
if pf && ~surprise
    time_index_constraint = maximum_lag + 1:maximum_lag + constrained_periods;
    second_system_size = pf * constrained_periods;
    J = zeros(second_system_size,second_system_size);
    r = zeros(second_system_size,1);
    indx_endo = zeros(second_system_size,1);
    col_count = 1;
    for j = 1:length(constrained_vars)
        indx_endo(col_count : col_count + constrained_periods - 1) = constrained_vars(j) + (time_index_constraint - 1) * ny;
        col_count = col_count + constrained_periods;
    end;
    make_ex_;
    make_y_;
    it = 1;
    convg = 0;
    while ~convg && it <= maxit
        disp('---------------------------------------------------------------------------------------------');
        disp(['iteration ' int2str(it)]);
        not_achieved = 1;
        alpha = 1;
        while not_achieved
            simul();
            result = sum(sum(isfinite(oo_.endo_simul(:,time_index_constraint)))) == ny * constrained_periods;
            if (~oo_.deterministic_simulation.status || ~result) && it > 1
                not_achieved = 1;
                alpha = alpha / 2;
                col_count = 1;
                for j = controlled_varexo
                    oo_.exo_simul(time_index_constraint,j) = (old_exo(:,j) + alpha * D_exo(col_count: col_count + constrained_periods - 1));
                    col_count = col_count + constrained_periods;
                end;
                disp(['Divergence in  Newton: reducing the path length alpha=' num2str(alpha)]);
                oo_.endo_simul = repmat(oo_.steady_state, 1, options_cond_fcst.periods + 2);
            else
                not_achieved = 0;
            end;
        end;
        
        y = oo_.endo_simul(constrained_vars(j), time_index_constraint);
        ys = oo_.endo_simul(indx_endo);
        col_count = 1;
        for j = 1:length(constrained_vars)
            y = oo_.endo_simul(constrained_vars(j), time_index_constraint);
            r(col_count:col_count+constrained_periods - 1) = (y - constrained_paths(j,1:constrained_periods))';
            col_count = col_count + constrained_periods;
        end;
        col_count = 1;
        for j = controlled_varexo
            for time = time_index_constraint
                saved = oo_.exo_simul(time,j);
                oo_.exo_simul(time,j) = oo_.exo_simul(time,j) + eps1;
                simul();
                J(:,col_count) = (oo_.endo_simul(indx_endo) - ys) / eps1;
                oo_.exo_simul(time,j) = saved;
                col_count = col_count + 1;
            end;
        end;
        normr = norm(r, 1);
        
        disp(['iteration ' int2str(it) ' error = ' num2str(normr)]);
        
        if normr <= eps
            convg = 1;
            disp('convergence achieved');
        else
            % Newton update on exogenous shocks
            old_exo = oo_.exo_simul(time_index_constraint,:);
            D_exo = - J \ r;
            col_count = 1;
            for j = controlled_varexo
                oo_.exo_simul(time_index_constraint,j) = (oo_.exo_simul(time_index_constraint,j) + D_exo(col_count: col_count + constrained_periods - 1));
                col_count = col_count + constrained_periods;
            end;
        end;
        it = it + 1;
    end;
    endo(maximum_lag + 1 :maximum_lag + options_cond_fcst.periods ,:) = oo_.endo_simul(:,maximum_lag + 1:maximum_lag + options_cond_fcst.periods )';
    exo = oo_.exo_simul;
else
    for t = 1:constrained_periods
        disp('=============================================================================================');
        disp(['t=' int2str(t) ' constrained_paths(:,t)=' num2str(constrained_paths(:,t)') ]);
        disp('=============================================================================================');
        if t == 1
            make_ex_;
            exo_init = oo_.exo_simul;
            make_y_;
        end;
        oo_.exo_simul = exo_init;
        oo_.endo_simul(:,1) = initial_conditions;
        convg = 0;
        it = 1;
        
        time_index_constraint = maximum_lag + 1:maximum_lag + constrained_periods - t + 1;
        
        second_system_size = surprise + pf * (constrained_periods - t + 1);
        indx_endo = zeros(second_system_size,1);
        col_count = 1;
        for j = 1:length(constrained_vars)
            if constrained_perfect_foresight(j)
                indx_endo(col_count : col_count + constrained_periods - t) = constrained_vars(j) + (time_index_constraint - 1) * ny;
                col_count = col_count + constrained_periods - t + 1;
            else
                indx_endo(col_count) = constrained_vars(j) + maximum_lag * ny;
                col_count = col_count + 1;
            end;
        end;
        
        J = zeros(second_system_size,second_system_size);
        
        r = zeros(second_system_size,1);
        
        while ~convg && it <= maxit
            disp('---------------------------------------------------------------------------------------------');
            disp(['iteration ' int2str(it) ' time ' int2str(t)]);
            not_achieved = 1;
            alpha = 1;
            while not_achieved
                simul();
                result = sum(sum(isfinite(oo_.endo_simul(:,time_index_constraint)))) == ny * length(time_index_constraint);
                if (~oo_.deterministic_simulation.status || ~result) && it > 1
                    not_achieved = 1;
                    alpha = alpha / 2;
                    col_count = 1;
                    for j = controlled_varexo
                        if constrained_perfect_foresight(j)
                            oo_.exo_simul(time_index_constraint,j) = (old_exo(time_index_constraint,j) + alpha * D_exo(col_count: col_count + constrained_periods - t));
                            col_count = col_count + constrained_periods - t + 1;
                        else
                            oo_.exo_simul(maximum_lag + 1,j) = old_exo(maximum_lag + 1,j) + alpha * D_exo(col_count);
                            col_count = col_count + 1;
                        end;
                    end;
                    disp(['Divergence in  Newton: reducing the path length alpha=' num2str(alpha) ' result=' num2str(result) ' sum(sum)=' num2str(sum(sum(isfinite(oo_.endo_simul(:,time_index_constraint))))) ' ny * length(time_index_constraint)=' num2str(ny * length(time_index_constraint)) ' oo_.deterministic_simulation.status=' num2str(oo_.deterministic_simulation.status)]);
                    oo_.endo_simul = [initial_conditions repmat(oo_.steady_state, 1, options_cond_fcst.periods + 1)];
                else
                    not_achieved = 0;
                end;
            end;
            if t==constrained_periods
                ycc = oo_.endo_simul;
            end;
            yc = oo_.endo_simul(:,maximum_lag + 1);
            ys = oo_.endo_simul(indx_endo);
            
            col_count = 1;
            for j = 1:length(constrained_vars)
                if constrained_perfect_foresight(j)
                    y = oo_.endo_simul(constrained_vars(j), time_index_constraint);
                    r(col_count:col_count+constrained_periods - t) = (y - constrained_paths(j,t:constrained_periods))';
                    col_count = col_count + constrained_periods - t + 1;
                else
                    y = yc(constrained_vars(j));
                    r(col_count) = y - constrained_paths(j,t);
                    col_count = col_count + 1;
                end;
            end
            
            disp('computation of derivatives w.r. to exogenous shocks');
            col_count = 1;
            for j = 1:length(controlled_varexo)
                j_pos = controlled_varexo(j);
                if constrained_perfect_foresight(j)
                    for time = time_index_constraint
                        saved = oo_.exo_simul(time,j_pos);
                        oo_.exo_simul(time,j_pos) = oo_.exo_simul(time,j_pos) + eps1;
                        simul();
                        J(:,col_count) = (oo_.endo_simul(indx_endo) - ys) / eps1;
                        oo_.exo_simul(time,j_pos) = saved;
                        col_count = col_count + 1;
                    end;
                else
                    saved = oo_.exo_simul(maximum_lag+1,j_pos);
                    oo_.exo_simul(maximum_lag+1,j_pos) = oo_.exo_simul(maximum_lag+1,j_pos) + eps1;
                    simul();
                    J(:,col_count) = (oo_.endo_simul(indx_endo) - ys) / eps1;
                    oo_.exo_simul(maximum_lag+1,j_pos) = saved;
                    col_count = col_count + 1;
                end;
            end;
            
            normr = norm(r, 1);
            
            disp(['iteration ' int2str(it) ' error = ' num2str(normr) ' at time ' int2str(t)]);
            
            if normr <= eps
                convg = 1;
                disp('convergence achieved');
            else
                % Newton update on exogenous shocks
                D_exo = - J \ r;
                old_exo = oo_.exo_simul;
                col_count = 1;
                for j = 1:length(controlled_varexo)
                    j_pos=controlled_varexo(j);
                    if constrained_perfect_foresight(j)
                        oo_.exo_simul(time_index_constraint,j_pos) = (oo_.exo_simul(time_index_constraint,j_pos) + D_exo(col_count: col_count + constrained_periods - t));
                        col_count = col_count + constrained_periods - t + 1;
                    else
                        oo_.exo_simul(maximum_lag + 1,j_pos) = oo_.exo_simul(maximum_lag + 1,j_pos) + D_exo(col_count);
                        col_count = col_count + 1;
                    end;
                end;
            end;
            it = it + 1;
        end;
        if ~convg
            error(['convergence not achived at time ' int2str(t) ' after ' int2str(it) ' iterations']);
        end;
        for j = 1:length(controlled_varexo)
            j_pos=controlled_varexo(j);
            if constrained_perfect_foresight(j)
                % in case of mixed surprise and perfect foresight
                % endogenous path at each date all the exogenous paths have to be
                % stored. The paths are stacked in exo.
                for time = time_index_constraint;
                    exo(past_val + time,j_pos) = oo_.exo_simul(time,j_pos);
                end
            else
                exo(maximum_lag + t,j_pos) = oo_.exo_simul(maximum_lag + 1,j_pos);
            end;
        end;
        past_val = past_val + length(time_index_constraint);
        if t < constrained_periods
            endo(maximum_lag + t,:) = yc;
        else
            endo(maximum_lag + t :maximum_lag + options_cond_fcst.periods ,:) = ycc(:,maximum_lag + 1:maximum_lag + options_cond_fcst.periods - constrained_periods + 1)';
        end;
        initial_conditions = yc;
    end;
end;
options_.periods = save_options_periods;
options_.dynatol.f = save_options_dynatol_f;
options_.initval_file = save_options_initval_file;
% disp('endo');
% disp(endo);
% disp('exo');
% disp(exo);

oo_.endo_simul = endo';
oo_.exo_simul = exo;
