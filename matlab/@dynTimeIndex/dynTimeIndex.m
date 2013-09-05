function t = dynTimeIndex() % --*-- Unitary tests --*--

% t = dynTimeIndex()
%
% Constructor for the dynTimeIndex class.
%
% INPUTS:
%  None.
%
% OUTPUTS:
%  * t, dynTimeIndex object.
%
% DESCRIPTION:
%  The dynTimeIndex object is used to shift backward or forward dynSeries objects. For instance, if ts
%  is a dynSeries object and t is a dynTimeIndex object then the following expressions are equivalent:
%
%      us = ts.lag()
%      us = ts.lag(1)
%      us = lag(ts,1)
%      us = ts(t-1)
%
%  This class has only one member: t = int8(0) when instantiated.

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

t = struct();

t.index = int8(0);

t = class(t,'dynTimeIndex');

%@test:1
%$ % Instantiate a dynTimeIndex object
%$ try
%$     u = dynTimeIndex();
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = isa(u,'dynTimeIndex');
%$ end
%$
%$ T = all(t);
%@eof:1