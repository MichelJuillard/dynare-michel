function ms_compute_probabilities(options_)
%function ms_simulation()
% calls ms sbvar probabilities mex file
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

compute_probabilities = create_probabilities_commandline(options_);
[err] = ms_sbvar_command_line(compute_probabilities);
mexErrCheck('ms_sbvar_command_line probabilities',err);

end

function opt=create_probabilities_commandline(options_)

opt = '-simulate';
opt = [opt set_flag_if_exists(options_.ms, 'output_file_tag')];
opt = [opt set_flag_if_exists(options_.ms, 'filtered_probabilities')];
opt = [opt set_flag_if_exists(options_.ms, 'output_file_tag')];

if options_.DynareRandomStreams.seed
    opt = [' -seed' num2str(options_.DynareRandomStreams.seed)];
end
end
