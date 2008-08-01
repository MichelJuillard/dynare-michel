function ldens = lpdfgbeta(x,a,b,aa,bb);

% function ldens = lpdfgbeta(x,a,b,aa,bb);
% log (generalized) BETA PDF 
%
% INPUTS
%    x:      density evatuated at x
%    a:      BETA distribution paramerer 
%    b:      BETA distribution paramerer 
%    aa:     lower bound 
%    bb:     upper bound 

% OUTPUTS
%    ldens:  the log (generalized) BETA PDF
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2008 Dynare Team
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

ldens = -betaln(a,b) + (a-1)*log(x-aa) + (b-1)*log(bb-x) - (a+b-1)*log(bb-aa);

%gammaln(a+b) - gammaln(a) - gammaln(b)
%betaln(a,b)
%pause
% 02/16/04 SA Interval [aa,bb] is the support of the probability density function. 
