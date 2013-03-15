function ps = pages(varargin)
%function ps = pages(varargin)
% Pages Class Constructor
%
% INPUTS
%   varargin        0 args  : empty pages object
%                   1 arg   : must be pages object (return a copy of arg)
%
% OUTPUTS
%   ps     [pages] pages object
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

switch nargin
    case 0
        ps = class(struct, 'pages', objArray());
    case 1
        assert(isa(varargin{1}, 'pages'), ...
               ['@pages.pages: With one arg to pages constructor, you must ' ...
                'pass an pages object or a char.']);
        ps = varargin{1};
    otherwise
        error('@pages.pages: invalid number of arguments');
end
end