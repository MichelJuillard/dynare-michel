function planner_objective_value = evaluate_planner_objective(M,options,oo)

%function oo1 = evaluate_planner_objective(dr,M,oo,options)
%  computes value of planner objective function     
% 
% INPUTS
%   M:        (structure) model description
%   options:  (structure) options
%   oo:       (structure) output results
%    
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2007-2012 Dynare Team
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

dr = oo.dr;
endo_nbr = M.endo_nbr;
exo_nbr = M.exo_nbr;
nstatic = dr.nstatic;
npred = dr.npred;
lead_lag_incidence = M.lead_lag_incidence;
beta = get_optimal_policy_discount_factor(M.params,M.param_names);
if options.ramsey_policy
    i_org = (1:M.orig_endo_nbr)';
else
    i_org = (1:M.endo_nbr)';
end
ipred = find(lead_lag_incidence(M.maximum_lag,:))';
order_var = dr.order_var;
LQ = true;

Gy = dr.ghx(nstatic+(1:npred),:);
Gu = dr.ghu(nstatic+(1:npred),:);
gy(dr.order_var,:) = dr.ghx;
gu(dr.order_var,:) = dr.ghu;

if options.ramsey_policy && options.order == 1 && ~options.linear
    options.order = 2;
    options.qz_criterium = 1+1e-6;
    [dr,info] = stochastic_solvers(oo.dr,0,M,options,oo);
    Gyy = dr.ghxx(nstatic+(1:npred),:);
    Guu = dr.ghuu(nstatic+(1:npred),:);
    Gyu = dr.ghxu(nstatic+(1:npred),:);
    Gss = dr.ghs2(nstatic+(1:npred),:);
    gyy(dr.order_var,:) = dr.ghxx;
    guu(dr.order_var,:) = dr.ghuu;
    gyu(dr.order_var,:) = dr.ghxu;
    gss(dr.order_var,:) = dr.ghs2;
    LQ = false;
end

ys = oo.dr.ys;

u = oo.exo_simul(1,:)';

[U,Uy,Uyy] = feval([M.fname '_objective_static'],ys,zeros(1,exo_nbr), ...
                   M.params);

Uyy = full(Uyy);

[Uyygygy, err] = A_times_B_kronecker_C(Uyy,gy,gy,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Uyygugu, err] = A_times_B_kronecker_C(Uyy,gu,gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Uyygygu, err] = A_times_B_kronecker_C(Uyy,gy,gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);

Wbar =U/(1-beta);
Wy = Uy*gy/(eye(npred)-beta*Gy);
Wu = Uy*gu+beta*Wy*Gu;
if LQ
    Wyy = Uyygygy/(eye(npred*npred)-beta*kron(Gy,Gy));
else
    Wyy = (Uy*gyy+Uyygygy+beta*Wy*Gyy)/(eye(npred*npred)-beta*kron(Gy,Gy));
end
[Wyygugu, err] = A_times_B_kronecker_C(Wyy,Gu,Gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wyygygu,err] = A_times_B_kronecker_C(Wyy,Gy,Gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
if LQ
    Wuu = Uyygugu+beta*Wyygugu;
    Wyu = Uyygygu+beta*Wyygygu;
    Wss = beta*Wuu*M.Sigma_e(:)/(1-beta);
else
    Wuu = Uy*guu+Uyygugu+beta*(Wy*Guu+Wyygugu);
    Wyu = Uy*gyu+Uyygygu+beta*(Wy*Gyu+Wyygygu);
    Wss = (Uy*gss+beta*(Wuu*M.Sigma_e(:)+Wy*Gss))/(1-beta);
end
% initialize yhat1 at the steady state
yhat1 = oo.steady_state;
if options.ramsey_policy
    % initialize le Lagrange multipliers to 0 in yhat2
    yhat2 = zeros(M.endo_nbr,1);
    yhat2(1:M.orig_endo_nbr) = oo.steady_state(1:M.orig_endo_nbr);
end
if ~isempty(M.endo_histval)
    % initialize endogenous state variable to histval if necessary
    yhat1(1:M.orig_endo_nbr) = M.endo_histval(1:M.orig_endo_nbr);
    if options.ramsey_policy
        yhat12(1:M.orig_endo_nbr) = M.endo_histval(1:M.orig_endo_nbr);
    end

end
yhat1 = yhat1(dr.order_var(nstatic+(1:npred)),1)-dr.ys(dr.order_var(nstatic+(1:npred)));
yhat2 = yhat2(dr.order_var(nstatic+(1:npred)),1)-dr.ys(dr.order_var(nstatic+(1:npred)));
u = oo.exo_simul(1,:)';

[Wyyyhatyhat1, err] = A_times_B_kronecker_C(Wyy,yhat1,yhat1,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wuuuu, err] = A_times_B_kronecker_C(Wuu,u,u,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wyuyhatu1, err] = A_times_B_kronecker_C(Wyu,yhat1,u,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
planner_objective_value(1) = Wbar+Wy*yhat1+Wu*u+Wyuyhatu1 ...
    + 0.5*(Wyyyhatyhat1 + Wuuuu+Wss);
if options.ramsey_policy
    [Wyyyhatyhat2, err] = A_times_B_kronecker_C(Wyy,yhat2,yhat2,options.threads.kronecker.A_times_B_kronecker_C);
    mexErrCheck('A_times_B_kronecker_C', err);
    [Wyuyhatu2, err] = A_times_B_kronecker_C(Wyu,yhat2,u,options.threads.kronecker.A_times_B_kronecker_C);
    mexErrCheck('A_times_B_kronecker_C', err);
    planner_objective_value(2) = Wbar+Wy*yhat2+Wu*u+Wyuyhatu2 ...
        + 0.5*(Wyyyhatyhat2 + Wuuuu+Wss);
end

if ~options.noprint
    disp(' ')
    disp('Approximated value of planner objective function')
    disp(['    - with initial Lagrange multipliers set to 0: ' ...
          num2str(planner_objective_value(2)) ])
    disp(['    - with initial Lagrange multipliers set to steady state: ' ...
          num2str(planner_objective_value(1)) ])
    disp(' ')
end