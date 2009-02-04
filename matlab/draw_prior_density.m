function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx,bayestopt_);
% function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx)
% Computes values of prior density at many points (before plotting)
%
% INPUTS
%    indx          [integer]    Parameter number.
%    bayestopt_    [structure]  Describes the prior beliefs.
%    
% OUTPUTS
%    x             [double]     Row vector, subset of 'abscissa' such as the density is less than 10
%    f             [double]     Row vector, subset of 'dens' such as the density is less than 10
%    abscissa      [double]     Row vector, abscissa 
%    dens          [double]     Row vector, density
%    binf:         [double]     Scalar, first element of x
%    bsup:         [double]     Scalar, last element of x
%    
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2004-2009 Dynare Team
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

pmean   = bayestopt_.pmean;
pshape  = bayestopt_.pshape; 
p1      = bayestopt_.p1;
p2      = bayestopt_.p2;
p3      = bayestopt_.p3;
p4      = bayestopt_.p4;

truncprior = 1e-3;
steps = 200;

switch pshape(indx)
 case 1  % Beta prior
    density = @(x,a,b,aa,bb) betapdf((x-aa)/(bb-aa), a, b)/(bb-aa);
    mu = (p1(indx)-p3(indx))/(p4(indx)-p3(indx));
    stdd = p2(indx)/(p4(indx)-p3(indx));
    a = (1-mu)*mu^2/stdd^2 - mu;
    b = a*(1/mu-1);
    aa = p3(indx);
    bb = p4(indx);
    infbound = betainv(truncprior,a,b)*(bb-aa)+aa;
    supbound = betainv(1-truncprior,a,b)*(bb-aa)+aa;
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b,aa,bb);
 case 2  % Generalized Gamma prior
    density = @(x,a,b,c) gampdf(x-c,a,b);
    mu = p1(indx)-p3(indx);
    b  = p2(indx)^2/mu;
    a  = mu/b;
    c = p3(indx);
    infbound = gaminv(truncprior,a,b)+c;
    supbound = gaminv(1-truncprior,a,b)+c;
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b,c);
 case 3  % Gaussian prior
    a = p1(indx);
    b = p2(indx);
    infbound = norminv(truncprior,a,b); 
    supbound = norminv(1-truncprior,a,b);
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = normpdf(abscissa,a,b);  
 case 4  % Inverse-gamma of type 1 prior
    nu = p2(indx);
    s  = p1(indx);
    infbound = 1/sqrt(gaminv(1-10*truncprior, nu/2, 2/s));
    supbound = 1/sqrt(gaminv(10*truncprior, nu/2, 2/s));
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = exp(lpdfig1(abscissa,s,nu));  
 case 5  % Uniform prior
    infbound = p1(indx); 
    supbound = p2(indx);
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = ones(1, steps) / (supbound-infbound);
 case 6  % Inverse-gamma of type 2 prior
    nu = p2(indx);
    s  = p1(indx);
    infbound = 1/(gaminv(1-10*truncprior, nu/2, 2/s));
    supbound = 1/(gaminv(10*truncprior, nu/2, 2/s));
    stepsize = (supbound-infbound)/steps;
    abscissa = infbound:stepsize:supbound;
    dens = exp(lpdfig2(abscissa,s,nu));
 otherwise
  error(sprintf('draw_prior_density: unknown distribution shape (index %d, type %d)', indx, pshape(indx)));
end 

k = [1:length(dens)];
if pshape(indx) ~= 5 
    [junk,k1] = max(dens);
    if k1 == 1 | k1 == length(dens)
        k = find(dens < 10);
    end
end
binf = abscissa(k(1));
bsup = abscissa(k(length(k)));
x = abscissa(k);
f = dens(k);