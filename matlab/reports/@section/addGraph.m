function o = addGraph(o, varargin)
%function o = addGraph(o, varargin)
% Add a graph to the Cell Array of graphs in the report
%
% INPUTS
%   1 args => add empty graph
%   2 args => add given graph
%   3 args => add graph at index
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

disp(['Processing Section Element: ' num2str(numElements(o)+1)]);
o.elements = o.elements.addGraph(varargin{:});
end
