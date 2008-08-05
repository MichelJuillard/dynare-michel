function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx);
% function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx)
% plots prior density
%
% INPUTS
%    indx:      parameter number
%    
% OUTPUTS
%    x:         subset of 'abscissa' such as the density is less than 10
%    f:         subset of 'dens' such as the density is less than 10
%    abscissa:  abscissa 
%    dens:      density
%    binf:      lower bound of the truncated prior
%    bsup:      upper bound of the truncated prior
%    
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2004-2008 Dynare Team
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

global bayestopt_

pmean   = bayestopt_.pmean;
pshape  = bayestopt_.pshape; 
p1      = bayestopt_.p1;
p2      = bayestopt_.p2;
p3      = bayestopt_.p3;
p4      = bayestopt_.p4;

truncprior = 10^(-3);

switch pshape(indx)
 case 1  % Beta prior
    density = inline('((bb-x).^(b-1)).*(x-aa).^(a-1)./(beta(a,b)*(bb-aa)^(a+b-1))','x','a','b','aa','bb');
    mu = (p1(indx)-p3(indx))/(p4(indx)-p3(indx));
    stdd = p2(indx)/(p4(indx)-p3(indx));
    a = (1-mu)*mu^2/stdd^2 - mu;
    b = a*(1/mu-1);
    aa = p3(indx);
    bb = p4(indx);
    infbound = betainv(truncprior,a,b)*(bb-aa)+aa;
    supbound = betainv(1-truncprior,a,b)*(bb-aa)+aa;
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b,aa,bb);
 case 2  % Generalized Gamma prior
    mu = p1(indx)-p3(indx);
    b  = p2(indx)^2/mu;
    a  = mu/b;
    infbound = gaminv(truncprior,a,b);
    supbound = gaminv(1-truncprior,a,b);
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = exp(lpdfgam(abscissa,a,b));
    abscissa = abscissa + p3(indx);
 case 3  % Gaussian prior
    density = inline('inv(sqrt(2*pi)*b)*exp(-0.5*((x-a)/b).^2)','x','a','b');
    a = p1(indx);
    b = p2(indx);
    infbound = norminv(truncprior,a,b); 
    supbound = norminv(1-truncprior,a,b);
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b);  
 case 4  % Inverse-gamma of type 1 prior
    density = inline('2*inv(gamma(nu/2))*(x.^(-nu-1))*((s/2)^(nu/2)).*exp(-s./(2*x.^2))','x','s','nu');
    nu = p2(indx);
    s  = p1(indx);
    a  = nu/2;
    b  = 2/s;
    infbound = 1/sqrt(gaminv(1-10*truncprior,a,b));
    supbound = 1/sqrt(gaminv(10*truncprior,a,b));
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,s,nu);  
 case 5  % Uniform prior
    density = inline('(x.^0)/(b-a)','x','a','b');
    a  = p1(indx);
    b  = p2(indx);
    infbound = a; 
    supbound = b;
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,a,b);  
 case 6  % Inverse-gamma of type 2 prior
    density = inline('inv(gamma(nu/2))*(x.^(-.5*(nu+2)))*((s/2)^(nu/2)).*exp(-s./(2*x))','x','s','nu');
    nu = p2(indx);
    s  = p1(indx);
    a  = nu/2;
    b  = 2/s;
    infbound = 1/(gaminv(1-truncprior,a,b));
    supbound = 1/(gaminv(truncprior,a,b));
    stepsize = (supbound-infbound)/200;
    abscissa = infbound:stepsize:supbound;
    dens = density(abscissa,s,nu);  
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