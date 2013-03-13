function o = addGraph(o, varargin)
%function o = addGraph(o, varargin)
% Add a section to the current page in the current report
%
% INPUTS
%   o        - report object
%   varargin - arguments to @section/addGraph.m
%
% OUTPUTS
%   o        - updated report object
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

lastPage = o.pages.numPages();
lastSection = o.pages(lastPage).numSections();
o.pages(lastPage).sections(lastSection) = ...
    o.pages(lastPage).sections(lastSection).addGraph(varargin{:});
end
