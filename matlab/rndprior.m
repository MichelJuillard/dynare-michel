function y = rndprior(bayestopt_)
% function y = rndprior(bayestopt_)
% Draws random number from the prior density
%
% INPUTS
%   bayestopt_:    structure characterizing priors
%    
% OUTPUTS
%   y:             drawn numbers vector              
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

pshape=bayestopt_.pshape;
pmean=bayestopt_.pmean;
p1=bayestopt_.p1;
p2=bayestopt_.p2;
p3=bayestopt_.p3;
p4=bayestopt_.p4;
 
for i=1:length(pmean),
    
    switch pshape(i)
        
     case 1 %'beta'
      mu = (pmean(i)-p3(i))/(p4(i)-p3(i));
      stdd = p2(i)/(p4(i)-p3(i));
      A = (1-mu)*mu^2/stdd^2 - mu;
      B = A*(1/mu - 1);
      y(1,i) = betarnd(A, B);
      y(1,i) = y(1,i) * (p4(i)-p3(i)) + p3(i);
      
     case 2 %'gamma'
      mu = pmean(i)-p3(i);
      B = p2(i)^2/mu;
      A = mu/B;
      y(1,i) = gamrnd(A, B) + p3(i);
      
     case 3 %'normal'
      MU = pmean(i);
      SIGMA = p2(i);
      y(1,i) = randn*SIGMA+ MU;
      
     case 4 %'invgamma'
      nu = p2(i);
      s = p1(i);
      y(1,i) = 1/sqrt(gamrnd(nu/2, 2/s));
      
     case 5 %'uniform'
      y(1,i) = rand*(p2(i)-p1(i)) + p1(i);
      
    end
end

% initial version by Marco Ratto

