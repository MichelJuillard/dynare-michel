function options_=set_ms_init_file(options_)
%function set_ms_init_file()
% Set options_.ms.init_file based on user input
%
% INPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
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

if ~isfield(options_.ms, 'initialization_file_tag')
    options_.ms.initialization_file_tag = options_.ms.output_file_tag;
end
options_.ms.init_file = ['init_' options_.ms.initialization_file_tag '.dat'];

if ~exist(options_.ms.init_file,'file')
    error(['ERROR: Could not find initialization file: ' options_.ms.init_file]);
end
end
