function o = series(varargin)
%function o = series(varargin)
% Series Class Constructor
%
% INPUTS
%   varargin        0 args  : empty series object
%                   1 arg   : must be series object (return a copy of arg)
%                   > 1 args: option/value pairs (see structure below for
%                   options)
%
% OUTPUTS
%   o   [series] series object
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

o = struct;

o.data = '';

o.graphLineColor = 'k';
o.graphLineStyle = '-';
o.graphLineWidth = 0.5;

o.graphMarker = '';
o.graphMarkerEdgeColor = 'auto';
o.graphMarkerFaceColor = 'auto';
o.graphMarkerSize = 6;

o.tableShowMarkers = false;
o.tableNegColor = 'red';
o.tablePosColor = 'blue';
o.tableMarkerLimit = 1e-4;

o.tableAlignRight = false;

if nargin == 1
    assert(isa(varargin{1}, 'series'),['@series.series: with one arg you ' ...
                        'must pass a series object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@series.series: options must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = fieldnames(o);

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        ind = strmatch(lower(pair{1}), lower(optNames), 'exact');
        assert(isempty(ind) || length(ind) == 1);
        if ~isempty(ind)
            o.(optNames{ind}) = pair{2};
        else
            error('@series.series: %s is not a recognized option.', pair{1});
        end
    end
end

% Create series object
o = class(o, 'series');
end