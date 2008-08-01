function lnprior = priordens(para, pshape, p1, p2, p3, p4)

% function lnprior = priordens(para, pshape, p1, p2, p3, p4)
% computes a prior density for the structural parameters of DSGE models
%
% INPUTS
%    para:       parameter value
%    pshape:     0 is point mass, both para and p2 are ignored
%                1 is BETA
%                2 is GAMMA
%                3 is NORMAL
%                4 is INVERTED GAMMA TYPE I
%                5 is UNIFORM
%                6 is INVERTED GAMMA TYPE II
%    p1:         mean
%    p2:         standard deviation
%    p3:         lower bound
%    p4:         upper bound
%    
% OUTPUTS
%    lnprior:    log of the prior density
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

lnprior = 0;
nprio 	= length(pshape);

i = 1;
while i <=  nprio;
    a = 0;	
    b = 0;
    if pshape(i) == 1;     % (generalized) BETA Prior 
        mu = (p1(i)-p3(i))/(p4(i)-p3(i));
        stdd = p2(i)/(p4(i)-p3(i));
        a = (1-mu)*mu^2/stdd^2 - mu;
        b = a*(1/mu - 1);
        lnprior = lnprior + lpdfgbeta(para(i),a,b,p3(i),p4(i))   ;
    elseif pshape(i) == 2; % GAMMA PRIOR 
     	b = p2(i)^2/(p1(i)-p3(i));
        a = (p1(i)-p3(i))/b;
        lnprior = lnprior + lpdfgam(para(i)-p3(i),a,b);
    elseif pshape(i) == 3; % GAUSSIAN PRIOR 
        lnprior = lnprior + lpdfnorm(para(i),p1(i),p2(i));
    elseif pshape(i) == 4; % INVGAMMA1 PRIOR 
     	lnprior = lnprior + lpdfig1(para(i),p1(i),p2(i));
    elseif pshape(i) == 5; % UNIFORM PRIOR 
     	lnprior = lnprior + log(1/(p2(i)-p1(i)));
    elseif pshape(i) == 6; % INVGAMMA2 PRIOR 
    	lnprior = lnprior + lpdfig2(para(i),p1(i),p2(i));
    end;
    i = i+1;
end;

% 10/11/03 MJ adapted from an earlier version in GAUSS by F. Schorfheide
%             and translated to Matlab by R. Wouters
% 11/18/03 MJ adopted M.Ratto's suggestion for inverse gamma
%             changed order of input parameters
% 01/16/04 MJ added invgamma2
%             for invgamma p2 is now standard error
% 16/02/04 SA changed beta prior call
