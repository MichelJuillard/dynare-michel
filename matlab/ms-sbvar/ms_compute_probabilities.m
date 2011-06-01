function [options_, oo_]=ms_compute_probabilities(M_, options_, oo_)
%function ms_simulation()
% Compute posterior mode regime probabilities
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
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

disp('Compute Marginal Data Density');
options_ = set_ms_estimation_flags_for_other_mex(options_);
options_ = set_ms_simulation_flags_for_other_mex(options_);
oo_ = set_oo_w_estimation_output(options_, oo_);

% setup command line options
opt = ['-probabilities -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -ft ' options_.ms.output_file_tag];

if options_.ms.filtered_probabilities
    opt = [opt ' -filtered' ];
    prob_out_file = ['filtered_' options_.ms.output_file_tag '.out'];
elseif options_.ms.real_time_smoothed_probabilities
    opt = [opt ' -real_time_smoothed' ];
    prob_out_file = 0;
else
    opt = [opt ' -smoothed' ];
    prob_out_file = ['smoothed_' options_.ms.output_file_tag '.out'];
end

% compute probabilities
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_sbvar_command_line probabilities',err);

% now we want to plot the probabilities for each chain
if ischar(prob_out_file)
    computed_probabilities = load(prob_out_file);
    plot_ms_probabilities(computed_probabilities,options_,M_.fname);
end
options_ = initialize_ms_sbvar_options(M_, options_);
end
