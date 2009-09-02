function resid(period)
% function resid(period)
% Computes residuals associated with the guess values
% 
% INPUTS
%    period
%    
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2009 Dynare Team
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

  global M_ options_ oo_ it_ z
  
  if M_.exo_nbr > 0
    oo_.exo_simul = ones(M_.maximum_lag+M_.maximum_lead+period,1)* ...
	oo_.exo_steady_state';
  end
  n = size(M_.lead_lag_incidence,2);
%  if ~ options_.initval_file | size(oo_.endo_simul,2) ~= period+M_.maximum_lag+M_.maximum_lead
  if ~ options_.initval_file 
    if size(oo_.steady_state,1) == 1 & oo_.steady_state == 0
      oo_.steady_state = zeros(size(oo_.steady_state,1),1) ;
    end
    oo_.endo_simul = oo_.steady_state*ones(1,period+M_.maximum_lag+M_.maximum_lead) ;
  end

  i = M_.lead_lag_incidence';
  iyr0 = find(i(:));

  y =oo_.endo_simul(:);
  z = zeros(n,period);
  fh = str2func([M_.fname '_dynamic']);
  for it_=M_.maximum_lag+1:period+M_.maximum_lag
        z(:,it_-M_.maximum_lag) = feval(fh,y(iyr0),oo_.exo_simul, M_.params, it_);
        iyr0 = iyr0 + n;
  end

  for i = 1:4
    disp(' ')
  end
  
  for i=1:length(z)
    if abs(z(i)) < options_.dynatol/100
        tmp = 0;
    else
        tmp = z(i);
    end
    disp(['Residual for equation number ' int2str(i) ' is equal to ' num2str(tmp)])
  end  
  for i = 1:2
    disp(' ')
  end