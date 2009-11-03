function info = perfect_foresight_simulation(init)
 % performs deterministic simulations with lead or lag on one period
 %
 % INPUTS
 %   none
 %
 % OUTPUTS
 %   none
 %
 % ALGORITHM
 %   Laffargue, Boucekkine, Juillard (LBJ)
 %   see Juillard (1996) Dynare: A program for the resolution and
 %   simulation of dynamic models with forward variables through the use
 %   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
 %
 % SPECIAL REQUIREMENTS
 %   None.
 
 % Copyright (C) 1996-2009 Dynare Team
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
 
 global M_ options_ oo_ 
 global ct_ it_ 
 
 persistent flag_init 
 persistent lead_lag_incidence dynamic_model ny nyp nyf nrs nrc iyf iyp isp is isf isf1 iz icf  
 
 if nargin==1 
     flag_init = [];
 end
 
 if isempty(flag_init) 
     lead_lag_incidence = M_.lead_lag_incidence; 
     dynamic_model = [M_.fname '_dynamic']; 
     ny   = size(oo_.endo_simul,1); 
     nyp  = nnz(lead_lag_incidence(1,:)); 
     nyf  = nnz(lead_lag_incidence(3,:)); 
     nrs  = ny+nyp+nyf+1; 
     nrc  = nyf+1; 
     iyf  = find(lead_lag_incidence(3,:)>0); 
     iyp  = find(lead_lag_incidence(1,:)>0); 
     isp  = 1:nyp; 
     is   = (nyp+1):(ny+nyp); 
     isf  = iyf+nyp; 
     isf1 = (nyp+ny+1):(nyf+nyp+ny+1);     
     iz   = 1:(ny+nyp+nyf); 
     icf  = 1:size(iyf,2); 
     flag_init = 1; 
     if nargin==1 
         return
     end
 end 
 
 endo_simul = oo_.endo_simul; 
 periods    = options_.periods; 
 
 stop    = 0 ; 
 it_init = M_.maximum_lag+1; 
 
 info.convergence = 1; 
 info.time  = 0; 
 info.error = 0; 
 info.iterations.time  = zeros(options_.maxit_,1); 
 info.iterations.error = info.iterations.time; 
 
 h1 = clock; 
 for iter = 1:options_.maxit_ 
     h2 = clock; 
     if ct_ == 0 
         c = zeros(ny*periods,nrc); 
     else
         c = zeros(ny*(periods+1),nrc);
     end
     it_ = it_init ; 
     z = [ endo_simul(iyp,it_-1) ; endo_simul(:,it_) ; endo_simul(iyf,it_+1) ]; 
     [d1,jacobian] = feval(dynamic_model,z,oo_.exo_simul, M_.params, it_); 
     jacobian = [jacobian(:,iz) , -d1]; 
     ic = 1:ny; 
     icp = iyp; 
     c(ic,:) = jacobian(:,is)\jacobian(:,isf1) ; 
     for it_ = it_init+(1:periods-1) 
         z = [ endo_simul(iyp,it_-1) ; endo_simul(:,it_) ; endo_simul(iyf,it_+1)]; 
         [d1,jacobian] = feval(dynamic_model,z,oo_.exo_simul, M_.params, it_); 
         jacobian = [jacobian(:,iz) , -d1]; 
         jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:); 
         ic = ic + ny; 
         icp = icp + ny; 
         c(ic,:) = jacobian(:,is)\jacobian(:,isf1); 
     end 
     if ct_ == 1 
         s = eye(ny);
         s(:,isf) = s(:,isf)+c(ic,1:nyf);
         ic = ic + ny;
         c(ic,nrc) = s\c(:,nrc);
         c = bksup0(c,ny,nrc,iyf,icf,periods);
         c = reshape(c,ny,periods+1);
         endo_simul(:,it_init+(0:periods)) = endo_simul(:,it_init+(0:periods))+options_.slowc*c;
     else 
         c = bksup0(c,ny,nrc,iyf,icf,periods); 
         c = reshape(c,ny,periods); 
         endo_simul(:,it_init+(0:periods-1)) = endo_simul(:,it_init+(0:periods-1))+options_.slowc*c; 
     end 
     err = max(max(abs(c./options_.scalv'))); 
     info.iterations.time(iter)  = etime(clock,h2); 
     info.iterations.error(iter) = err;   
     if err < options_.dynatol 
         stop = 1; 
         info.time  = etime(clock,h1); 
         info.error = err; 
         info.iterations.time = info.iterations.time(1:iter); 
         info.iterations.error  = info.iterations.error(1:iter); 
         oo_.endo_simul = endo_simul; 
         break 
     end
 end 
 
 if ~stop 
     info.time  = etime(clock,h1);
     info.error = err;
     info.convergence = 0;
 end