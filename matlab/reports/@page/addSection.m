function p = addSection(p, varargin)
%function p = addSection(p, varargin)
% Add a section to the Cell Array of sections in the report
%
% INPUTS
%   1 args => add empty section
%   2 args => add given section
%   3 args => add section at index
%
% OUTPUTS
%   updated report object
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

assert(nargin >= 1 && nargin <= 3, ['incorrect number of arguments passed ' ...
    'to addSection']);
assert(isa(p, 'page'), 'First argument must be a page object');
if nargin > 1
    assert(isa(varargin{1},'section'), ['Optional 2nd arg to addSection must be a ' ...
        'Section']);
    if nargin > 2
        assert(isnumeric(varargin{2}), ['Optional 3rd arg to addSection must be ' ...
            'an index']);
    end
end

if nargin == 1
    p.sections = p.sections.addSection(section());
elseif nargin == 2 || nargin == 3
    p.sections = p.sections.addSection(varargin{:});
end
end
