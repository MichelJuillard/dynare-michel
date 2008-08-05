function bounds = prior_bounds(bayestopt)
% function bounds = prior_bounds(bayestopt)
% computes practical bounds for prior density
%
% INPUTS
%    bayestopt:        structure characterizing priors (shape, mean, p1..p4)
%    
% OUTPUTS
%    bounds:           matrix specifying bounds (row= parameter, column=upper&lower bound)
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

global options_

pshape = bayestopt.pshape;
pmean = bayestopt.pmean;
p1 = bayestopt.p1;
p2 = bayestopt.p2;
p3 = bayestopt.p3;
p4 = bayestopt.p4;

n = length(pmean);
bounds = zeros(n,2);

for i=1:n
  switch pshape(i)
    case 1
      mu = (pmean(i)-p3(i))/(p4(i)-p3(i));
      stdd = p2(i)/(p4(i)-p3(i));
      A = (1-mu)*mu^2/stdd^2 - mu;
      B = A*(1/mu - 1);
      bounds(i,1) = betainv(options_.prior_trunc,A,B)*(p4(i)-p3(i))+p3(i);
      bounds(i,2) = betainv(1-options_.prior_trunc,A,B)*(p4(i)-p3(i))+p3(i);
    case 2
      b = p2(i)^2/(pmean(i)-p3(i));
      a = (pmean(i)-p3(i))/b;
      bounds(i,1) = gaminv(options_.prior_trunc,a,b)+p3(i);
      bounds(i,2) = gaminv(1-options_.prior_trunc,a,b)+p3(i);
    case 3
      bounds(i,1) = norminv(options_.prior_trunc,pmean(i),p2(i));
      bounds(i,2) = norminv(1-options_.prior_trunc,pmean(i),p2(i));
    case 4
      bounds(i,1) = 1/sqrt(gaminv(1-options_.prior_trunc, p2(i)/2, 2/p1(i)));
      bounds(i,2) = 1/sqrt(gaminv(options_.prior_trunc, p2(i)/2, 2/p1(i)));
    case 5
      bounds(i,1) = p1(i);
      bounds(i,2) = p2(i);
    case 6
      bounds(i,1) = 1/gaminv(1-options_.prior_trunc, p2(i)/2, 2/p1(i));
      bounds(i,2) = 1/gaminv(options_.prior_trunc, p2(i)/2, 2/p1(i));
    otherwise
      error(sprintf('prior_bounds: unknown distribution shape (index %d, type %d)', i, pshape(i)));
  end
end
