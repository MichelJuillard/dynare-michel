function [info,endo_simul] = perfect_foresight_simulation(endo_simul,exo_simul,compute_linear_solution,steady_state)
% Performs deterministic simulations with lead or lag on one period
%
% INPUTS
%   endo_simul                  [double]     n*T matrix, where n is the number of endogenous variables.
%   exo_simul                   [double]     q*T matrix, where q is the number of shocks.
%   compute_linear_solution     [integer]    scalar equal to zero or one.     
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

global M_ options_ it_

persistent flag_init 
persistent lead_lag_incidence dynamic_model ny nyp nyf nrs nrc iyf iyp isp is isf isf1 iz icf  
persistent ghx

if isempty(flag_init) 
    lead_lag_incidence = M_.lead_lag_incidence; 
    dynamic_model = [M_.fname '_dynamic'];
    ny   = size(endo_simul,1); 
    nyp  = nnz(lead_lag_incidence(1,:));% number of lagged variables.
    nyf  = nnz(lead_lag_incidence(3,:));% number of leaded variables. 
    nrs  = ny+nyp+nyf+1; 
    nrc  = nyf+1; 
    iyf  = find(lead_lag_incidence(3,:)>0);% indices for leaded variables. 
    iyp  = find(lead_lag_incidence(1,:)>0);% indices for lagged variables. 
    isp  = 1:nyp;
    is   = (nyp+1):(nyp+ny); % Indices for contemporaneaous variables. 
    isf  = iyf+nyp; 
    isf1 = (nyp+ny+1):(nyf+nyp+ny+1);     
    iz   = 1:(ny+nyp+nyf);
    icf  = 1:size(iyf,2);
    flag_init = 1;
end

if nargin<3
    compute_linear_solution = 0;
    if nargin<4
        error('The steady state (fourth input argument) is missing!');
    end
end

if compute_linear_solution
    [dr,info]=resol(steady_state,0);
    ghx = dr.ghx(end-dr.nfwrd+1:end,:);
end

periods = options_.periods; 

stop    = 0 ; 
it_init = M_.maximum_lag+1; 

info.convergence = 1; 
info.time  = 0; 
info.error = 0; 
info.iterations.time  = zeros(options_.maxit_,1); 
info.iterations.error = info.iterations.time; 

last_line = options_.maxit_;
error_growth = 0;

h1 = clock;
for iter = 1:options_.maxit_ 
    h2 = clock;
    if options_.terminal_condition
        c = zeros(ny*(periods+1),nrc);
    else
        c = zeros(ny*periods,nrc);
    end
    it_ = it_init;
    z = [ endo_simul(iyp,it_-1) ; endo_simul(:,it_) ; endo_simul(iyf,it_+1) ]; 
    [d1,jacobian] = feval(dynamic_model,z,exo_simul, M_.params, it_); 
    jacobian = [jacobian(:,iz) , -d1]; 
    ic = 1:ny;
    icp = iyp;
    c(ic,:) = jacobian(:,is)\jacobian(:,isf1) ; 
    for it_ = it_init+(1:periods-1)
        z = [ endo_simul(iyp,it_-1) ; endo_simul(:,it_) ; endo_simul(iyf,it_+1)]; 
        [d1,jacobian] = feval(dynamic_model,z,exo_simul, M_.params, it_); 
        jacobian = [jacobian(:,iz) , -d1]; 
        jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:); 
        ic = ic + ny;
        icp = icp + ny;
        c(ic,:) = jacobian(:,is)\jacobian(:,isf1); 
    end
    if options_.terminal_condition
        if options_.terminal_condition==1% Terminal condition is Y_{T} = Y_{T+1} 
            s = eye(ny);
            s(:,isf) = s(:,isf)+c(ic,1:nyf);
            ic = ic + ny;
            c(ic,nrc) = s\c(ic,nrc);
        else% Terminal condition is Y_{T}-Y^{\star} = TransitionMatrix*(Y_{T+1}-Y^{\star})
            z = [ endo_simul(iyp,it_-1) ; endo_simul(:,it_) ; endo_simul(iyf,it_+1) ] ;
            [d1,jacobian] = feval(dynamic_model,z,exo_simul, M_.params, it_);
            jacobian = [jacobian(:,iz) -d1];
            jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:) ;
            ic = ic + ny;
            icp = icp + ny;
            s = jacobian(:,is);
            s(:,iyp-nyp) = s(:,iyp-nyp)+jacobian(:,isf)*ghx;
            c (ic,:) = s\jacobian(:,isf1);
        end
        c = bksup0(c,ny,nrc,iyf,icf,periods);
        c = reshape(c,ny,periods+1);
        endo_simul(:,it_init+(0:periods)) = endo_simul(:,it_init+(0:periods))+options_.slowc*c;
    else% Terminal condition is Y_{T}=Y^{\star}
        c = bksup0(c,ny,nrc,iyf,icf,periods);
        c = reshape(c,ny,periods);
        endo_simul(:,it_init+(0:periods-1)) = endo_simul(:,it_init+(0:periods-1))+options_.slowc*c; 
    end
    err = max(max(abs(c))); 
    info.iterations.time(iter)  = etime(clock,h2); 
    info.iterations.error(iter) = err;
    % if iter>1
    %     error_growth = error_growth + (info.iterations.error(iter)>info.iterations.error(iter-1));
    % end
    % if isnan(err) || error_growth>3
    %     last_line = iter;
    %     break
    % end
    if err < options_.dynatol
        stop = 1;
        info.time  = etime(clock,h1); 
        info.error = err;
        info.iterations.time = info.iterations.time(1:iter); 
        info.iterations.error  = info.iterations.error(1:iter);
        break
    end
end

if options_.terminal_condition==2
    distance_to_steady_state = abs(((endo_simul(:,end-1)-endo_simul(:,end))./endo_simul(:,end)))*100;
    disp('Distance to steady state at the end is (in percentage):')
    distance_to_steady_state
end

if ~stop
    info.time  = etime(clock,h1);
    info.error = err;
    info.convergence = 0;
    info.iterations.time  = info.iterations.time(1:last_line);
    info.iterations.error = info.iterations.error(1:last_line);
    info.iterations.error
    endo_simul = [ ];
end