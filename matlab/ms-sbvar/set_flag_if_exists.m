function flagstr = set_flag_if_exists(msoptions, flag)
%function set_flag_if_exists()
% adds a flag to the command line if it exists in the msoptions structure
%
% INPUTS
%    msoptions:  options_.ms structure
%    flag:       string representing the flag to check
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011 Dynare Team
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

if ~isfield(msoptions, flag)
    flagstr = '';
    return;
end

% Convert from Dynare options to ms-sbvar c-code options here
% to allow for clear programming in the rest of the Dynare
% matlab files
switch flag
    case 'output_file_tag'              % general option
        cFlag = 'ft';
    case 'convergence_starting_value'   % estimation options
        cFlag = 'cb';
    case 'convergence_ending_value'
        cFlag = 'ce';
    case 'convergence_increment_value'
        cFlag = 'ci';
    case 'max_iterations_starting_value'
        cFlag = 'ib';
    case 'max_iterations_increment_value'
        cFlag = 'ii';
    case 'max_block_iterations'
        cFlag = 'mb';
    case 'max_repeated_optimization_runs'
        cFlag = 'repeat_max';
    case 'function_convergence_criterion'
        cFlag = 'repeat_tol_obj';
    case 'parameter_convergence_criterion'
        cFlag = 'repeat_tol_params';
    case 'number_of_large_perturbations'
        cFlag = 'random';
    case 'number_of_small_perturbations'
        cFlag = 'random_small';
    case 'number_of_posterior_draws_after_perturbation'
        cFlag = 'random_small_ndraws';
    case 'max_number_of_stages'
        cFlag = 'random_max';
    case 'random_function_convergence_criterion'
        cFlag = 'random_tol_obj';
    case 'random_parameter_convergence_criterion'
        cFlag = 'random_tol_params';
    case 'mh_replic'                    % simulation options
        cFlag = 'ndraws';
    case 'drop'
        cFlag = 'burnin';
    case 'thinning_factor'
        cFlag = 'thin';
    case 'adaptive_mh_draws'
        cFlag = 'mh';
    case 'load_mh_file'                 % mdd options
        cFlag = 'pf';
    case 'mdd_proposal_draws'
        cFlag = 'mdd_proposal_draws';
    case 'mdd_use_mean_center'
        cFlag = 'use_mean';
    case 'filtered_probabilities'       % probabilities options
        cFlag = 'filtered';
    case 'real_time_smoothed_probabilities'
        cFlag = 'real_time_smoothed';
    otherwise
        cFlag = flag;
end

flagstr = [' -' cFlag ' ' num2str(eval(['msoptions.' flag]))];

end