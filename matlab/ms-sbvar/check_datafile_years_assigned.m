function check_datafile_years_assigned(options_)
% function check_datafile_years_assigned(options_)
% check that datafile, initial_year and final_year were assigned
%
% INPUTS
%    options_:    (struct)    options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2012 Dynare Team
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

if isempty(options_.ms.initial_year)
    error('Must set initial_year option');
end

if isempty(options_.ms.final_year)
    error('Must set final_year option');
end

if isempty(options_.datafile)
    error('Must set datafile option');
end
