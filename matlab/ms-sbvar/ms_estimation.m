function ms_estimation(options_)
%function ms_estimation()
% calls ms sbvar estimation mex file
%
% INPUTS
%    options_
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

% cleanup old files
clean_ms_files(options_.ms.output_file_tag);

% general setup
ms_sbvar_setup(options_);

% estimation
perform_estimation = create_estimation_commandline(options_);
[err] = ms_sbvar_command_line(perform_estimation);
mexErrCheck('ms_sbvar_command_line estimation',err);

end

function opt=create_estimation_commandline(options_)

opt = '-estimate';
opt = [opt set_flag_if_exists(options_.ms, 'output_file_tag')];
opt = [opt set_flag_if_exists(options_.ms, 'convergence_starting_value')];
opt = [opt set_flag_if_exists(options_.ms, 'convergence_ending_value')];
opt = [opt set_flag_if_exists(options_.ms, 'convergence_increment_value')];
opt = [opt set_flag_if_exists(options_.ms, 'max_iterations_starting_value')];
opt = [opt set_flag_if_exists(options_.ms, 'max_iterations_increment_value')];
opt = [opt set_flag_if_exists(options_.ms, 'max_block_iterations')];
opt = [opt set_flag_if_exists(options_.ms, 'max_repeated_optimization_runs')];
opt = [opt set_flag_if_exists(options_.ms, 'function_convergence_criterion')];
opt = [opt set_flag_if_exists(options_.ms, 'parameter_convergence_criterion')];
opt = [opt set_flag_if_exists(options_.ms, 'number_of_large_perturbations')];
opt = [opt set_flag_if_exists(options_.ms, 'number_of_small_perturbations')];
opt = [opt set_flag_if_exists(options_.ms, 'number_of_posterior_draws_after_perturbation')];
opt = [opt set_flag_if_exists(options_.ms, 'max_number_of_stages')];
opt = [opt set_flag_if_exists(options_.ms, 'random_function_convergence_criterion')];
opt = [opt set_flag_if_exists(options_.ms, 'random_parameter_convergence_criterion')];

if options_.DynareRandomStreams.seed
    opt = [' -seed' num2str(options_.DynareRandomStreams.seed)];
end

% if ischar(options_.ms.mode_file) > 0 && exist(options_.ms.mode_file) > 0
%     opt = [opt ' -fp ' options_.ms.mode_file];
% end

end
