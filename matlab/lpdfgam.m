function  ldens = lpdfgam(x,a,b);

% function ldens = lpdfgam(x,a,b)
% log GAMMA PDF
%
% INPUTS
%    x:      density evatuated at x
%    a:      GAMMA distribution paramerer 
%    b:      GAMMA distribution paramerer 
%    
% OUTPUTS
%    ldens:  the log GAMMA PDF
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

  ldens = -gammaln(a) -a*log(b)+ (a-1)*log(x) -x/b ;

% 10/11/03  MJ adapted from an earlier GAUSS version by F. Schorfeide,
%              translated to MATLAB by R. Wouters  
%              use MATLAB gammaln rather than lngam
