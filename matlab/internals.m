function internals(flag, varargin)

%@info:
%! @deftypefn {Function File} internals (@var{flag},@var{a},@var{b}, ...)
%! @anchor{internals}
%! @sp 1
%! This command provides internal documentation and unitary tests for the matlab routines.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item flag
%! Mandatory argument: --doc (for displaying internal documentation) or --test (for performing unitary tests).
%! @item b
%! Name of the routine to be tested or for which internal documentation is needed.
%! @item c
%! Name of the routine to be tested.
%! @item d
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! None.
%! @sp 2
%! @strong{Examples}
%! @sp 1
%! The following instruction:
%! @sp 1
%! @example
%! internals --info particle/local_state_iteration
%! @end example
%! will display the internal documentation of the routine local_state_iteration located in the particle subfolder of the matlab directory.
%! @sp 1
%! The following instruction:
%! @sp 1
%! @example
%! internals --test particle/local_state_iteration
%! @end example
%! will execute the unitary tests associated the routine local_state_iteration.
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! [1] It is not possible to display the internal documentation of more than one routine.
%! @sp 1
%! [2] It is possible to perform unitary tests on a list of routines.
%! @sp 1
%! [3] For displaying the internal documentation, matlab calls texinfo which has to be installed.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! None.
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{utilities/tests/dynTest} @ref{utilities/doc/dynInfo}
%! @end deftypefn
%@eod:

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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

more off

if strcmpi(flag,'--test')
    if nargin>1
        dynare_config([],0);
        number_of_matlab_routines = length(varargin);
        for i=1:number_of_matlab_routines
            dynTest(varargin{i});
        end
    else
        disp('You have to specify at least one matlab routine after --test flag!')
    end
    return
end

if strcmpi(flag,'--info')
    if nargin==2
        dynare_config([],0);
        dynInfo(varargin{1})
    else
        if nargin<2
            disp('You have to specify a matlab routine after --info flag!')
        else
            disp('I can only show internal documentation for one matlab routine!')
        end
    end
    return
end

disp('You should read the manual...')
