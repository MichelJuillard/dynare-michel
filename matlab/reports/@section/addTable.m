function o = addTable(o, varargin)
%function o = addTable(o, varargin)
% Add a table to the Cell Array of tables in the report
%
% INPUTS
%   1 args => add empty table
%   2 args => add given table
%   3 args => add table at index
%
% OUTPUTS
%   updated section object
%
% SPECIAL REQUIREMENTS
%   none

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

assert(isa(o, 'section'), 'First argument must be a section object');
if nargin == 1
    o.elements = o.elements.addTable();
elseif nargin == 2 || nargin == 3
    o.elements = o.elements.addTable(varargin{:});
end
end
