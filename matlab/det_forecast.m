function det_cond_forecast(constrained_paths, constrained_vars, options_cond_fcst, constrained_perfect_foresight)
% Computes forecasts using the schocks retrieved from a condition forecast for a deterministic model.
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
    options_cond_fcst.periods = 60;
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
        if strcmp(exo_names(i,:), options_cond_fcst.controlled_varexo(j,:))
            controlled_varexo(j) = i;
        end
    end
end

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
maxit = options_.simul.maxit;

% Check the solution using a unconditional forecast (soft tune)

initial_conditions = oo_.steady_state;
terminal_conditions = oo_.steady_state;
exo = oo_.exo_simul;
T = options_.periods + 2;
endo_simul = zeros(ny, T);
endo_simul(:,1) = initial_conditions;
endo_simul(:,T) = initial_conditions;
exo_simul = zeros(T, nx);
exo_simul(1,:) = [oo_.exo_steady_state' oo_.exo_det_steady_state'];
exo_simul(T,:) = [oo_.exo_steady_state'  oo_.exo_det_steady_state'];
past_val = 0;

if pf && ~surprise
    make_ex_;
    make_y_;
    oo_.endo_simul(:,1) = initial_conditions;
    oo_.endo_simul(:,options_.periods + 2) = terminal_conditions;
    %oo_.exo_simul = repmat(oo_.exo_steady_state, options_.periods + 2, 1);
    oo_.exo_simul = exo;
    simul();
    endo_simul = oo_.endo_simul;
    exo_simul = oo_.exo_simul;
else
    for t=1:constrained_periods
        make_ex_;
        make_y_;
        disp(['t=' int2str(t) ' constrained_periods=' int2str(constrained_periods)]);
        oo_.endo_simul(:,1) = initial_conditions;
        oo_.endo_simul(:,options_.periods + 2) = terminal_conditions;
        time_index_constraint = maximum_lag + 1:maximum_lag + constrained_periods - t + 1;
        if t <= constrained_periods
            for j = controlled_varexo
                if constrained_perfect_foresight(j)
                    for time = time_index_constraint;
                        oo_.exo_simul(time,j) = exo(past_val + time,j);
                    end;
                    oo_.exo_simul(time+1, j)= oo_.exo_steady_state(j);
                else
                    oo_.exo_simul(maximum_lag + 1,j) = exo(maximum_lag + t,j);
                end;
            end;
        else
            tmp = size(oo_.exo_simul,1);
            oo_.exo_simul = repmat(oo_.exo_steady_state',tmp,1);
        end;
        past_val = past_val + length(time_index_constraint);
        simul();
        initial_conditions = oo_.endo_simul(:,2);
        if t < constrained_periods
            endo_simul(:,t+1) = initial_conditions;
            exo_simul(t+1,:) = oo_.exo_simul(2,:);
        else
            endo_simul(:,t + 1:t + options_cond_fcst.periods + maximum_lead) = oo_.endo_simul(:,maximum_lag + 1:maximum_lag + options_cond_fcst.periods + maximum_lead);
            exo_simul(t+1,:) = oo_.exo_simul(2,:);
        end;
    end;
end;
oo_.endo_simul = endo_simul;
oo_.exo_simul = exo_simul;
