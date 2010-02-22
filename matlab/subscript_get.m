function deriv = subscript_get(nargsout, func, args, varargin)
% function deriv = subscript_get(nargsout, func, args, varargin)
% returns the appropriate entry of the return argument from a  user-defined function
% which returns either the jacobian or hessian (or both)
%
% INPUTS
%    nargsout   [int]              integer indicating the number of the return argument containing the jac/hess
%    func       [function handle]  associated with the function
%    args       [cell array]       arguments provided to func
%    varargin   [cell array]       arguments showing the index (or indices) of the element to be returned
%
% OUTPUTS
%    deriv      [double]           the (element1,element2) entry of the hessian
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010 Dynare Team
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

switch size(varargin,2)
    case 1 %first deriv
        switch nargsout
            case 1
                [outargs{1}] = func(args{:});
                deriv = outargs{1}(varargin{1});
            case 2
                [outargs{1} outargs{2}] = func(args{:});
                deriv = outargs{2}(varargin{1});
            otherwise
                error('Wrong number of output arguments (%d) passed to subscript_get().',nargsout);
        end
    case 2 %second deriv
        switch nargsout
            case 1
                [outargs{1}] = func(args{:});
                deriv = outargs{1}(varargin{1},varargin{2});
            case 3
                [outargs{1} outargs{2} outargs{3}] = func(args{:});
                deriv = outargs{3}(varargin{1},varargin{2});
            otherwise
                error('Wrong number of output arguments (%d) passed to subscript_get().',nargsout);
        end
    otherwise
        error('Wrong number of indices (%d) was passed to subscript_get().',size(varargin,2));
end
end
