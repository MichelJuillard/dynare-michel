function options_=set_input_options(options_)
%function set_input_options()
% get mode_file from filename if it's not in mode_file
%
% INPUTS
%    options_
%
% OUTPUTS
%    options_
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

% init file
if ~isfield(options_.ms, 'init_file')
    options_.ms.init_file = ['init_' options_.ms.output_file_tag '.dat'];
end

if ~exist(options_.ms.init_file, 'file')
    error(['ERROR: cannot find init_file: ' options_.ms.init_file]);
end

% free parameter file
if ~isfield(options_.ms, 'free_param_file')
    options_.ms.free_param_file = ['est_free_' options_.ms.output_file_tag '.out'];
end

if ~exist(options_.ms.free_param_file, 'file')
    error(['ERROR: cannot find free_param_file: ' options_.ms.free_param_file]);
end

% simulation file
if ~isfield(options_.ms, 'load_mh_file')
    options_.ms.load_mh_file = ['simulation_' options_.ms.output_file_tag '.out'];
end

if ~exist(options_.ms.load_mh_file, 'file')
    error(['ERROR: Could not find Metropolis Hastings file: ' options_.ms.load_mh_file]);
end
end
