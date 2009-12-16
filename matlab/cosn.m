function co = cosn(H);

% function co = cosn(H);
% computes the cosine of the angle between the H(:,1) and its
% projection onto the span of H(:,2:end)
%
% Not the same as multiple correlation coefficient since the means are not
% zero
%
% Copyright (C) 2008 Dynare Team
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

y = H(:,1);
X = H(:,2:end);

% y = H(:,1);
% X = H(:,2:end);

yhat =  X*(X\y);
if rank(yhat),
    co = y'*yhat/sqrt((y'*y)*(yhat'*yhat));
else
    co=0;
end



