function dyn_saveas(h,fname,DynareOptions)
%function dyn_saveas(h,fname,DynareOptions)
% save figures for DYNARE
%
% INPUTS
%    h     : figure handle
%    fname : name of the saved figure
%    DynareOptions: dynare options
%
% OUTPUTS
%    none
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

if (strcmp('eps',DynareOptions.graph_format))
    eval(['print -depsc2 ' fname '.eps']);
end
if (strcmp('pdf',DynareOptions.graph_format))
    eval(['print -dpdf ' fname]);
end
if ~exist('OCTAVE_VERSION') ...
   && (strcmp('fig',DynareOptions.graph_format))
    if DynareOptions.nodisplay
        set(h, 'Visible','on');
    end
    saveas(h,[fname '.fig']);
end
if DynareOptions.nodisplay
    close(h);
end
