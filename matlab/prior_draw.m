function pdraw = prior_draw(init,cc)
% function pdraw = prior_draw(init,cc)
% Builds one draw from the prior distribution. 
% 
% INPUTS
%   o init           [integer]  scalar equal to 1 (first call) or 0.
%   o cc             [double]   two columns matrix (same as in
%                               metropolis.m), constraints over the
%                               parameter space (upper and lower bounds).
%    
% OUTPUTS 
%   o pdraw          [double]   draw from the joint prior density. 
%
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2006-2008 Dynare Team
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

global estim_params_  bayestopt_
persistent fname npar bounds pshape pmean pstd a b p1 p2 p3 p4 condition
  
if init
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    pshape = bayestopt_.pshape;
    pmean = bayestopt_.pmean;
    pstd  = bayestopt_.pstdev;
    p1 = bayestopt_.p1;
    p2 = bayestopt_.p2;
    p3 = bayestopt_.p3;
    p4 = bayestopt_.p4;
    a = zeros(npar,1);
    b = zeros(npar,1); 
    if nargin == 2
        bounds = cc;
    else
        bounds = kron(ones(npar,1),[-Inf Inf]);
    end
    for i = 1:npar
        switch pshape(i)
          case 3% Gaussian prior
            b(i) = pstd(i)^2/(pmean(i)-p3(i));
            a(i) = (pmean(i)-p3(i))/b(i);
          case 1% Beta prior
            mu = (p1(i)-p3(i))/(p4(i)-p3(i));
            stdd = p2(i)/(p4(i)-p3(i));
            a(i) = (1-mu)*mu^2/stdd^2 - mu;
            b(i) = a(i)*(1/mu - 1);        
          case 2;%Gamma prior
            mu = p1(i)-p3(i);
            b(i) = p2(i)^2/mu;
            a(i) = mu/b(i);
          case {5,4,6}
            % Nothing to do here
            %
            % 4: Inverse gamma, type 1, prior
            %    p2(i) = nu
            %    p1(i) = s
            % 6: Inverse gamma, type 2, prior
            %    p2(i) = nu
            %    p1(i) = s
            % 5: Uniform prior
            %    p3(i) and p4(i) are used.
          otherwise
           error(sprintf('prior_draw: unknown distribution shape (index %d, type %d)', i, pshape(i)));
        end
        pdraw = zeros(npar,1);  
    end
    condition = 1;
    pdraw = zeros(npar,1);
    return
end


for i = 1:npar
    switch pshape(i)
      case 5% Uniform prior.
        pdraw(i) = rand*(p4(i)-p3(i)) + p3(i);
      case 3% Gaussian prior.
        while condition
            tmp = randn*pstd(i) + pmean(i);
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end
      case 2% Gamma prior.
        while condition
            g = gamrnd(a(i),b(i)) + p3(i);
            if g >= bounds(i,1) && g <= bounds(i,2)
                pdraw(i) = g;
                break
            end
        end
      case 1% Beta distribution
        while condition
            tmp = betarnd(a(i), b(i));
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = p3(i)+tmp*(p4(i)-p3(i));
                break
            end
        end
      case 4% INV-GAMMA1 distribution
        while condition
            tmp = sqrt(1/gamrnd(p2(i)/2,2/p1(i)));
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end        
      case 6% INV-GAMMA2 distribution  
        while condition
            tmp = 1/gamrnd(p2(i)/2,2/p1(i));
            if tmp >= bounds(i,1) && tmp <= bounds(i,2)
                pdraw(i) = tmp;
                break
            end
        end
      otherwise
       error(sprintf('prior_draw: unknown distribution shape (index %d, type %d)', i, pshape(i)));
    end
end
