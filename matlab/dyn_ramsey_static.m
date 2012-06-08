function [steady_state,params,check] = dyn_ramsey_static(x,M,options_,oo)

% function  [steady_state,params,check] = dyn_ramsey_static_(x)
% Computes the static first order conditions for optimal policy
%
% INPUTS
%    x:         vector of endogenous variables or instruments
%
% OUTPUTS
%    resids:    residuals of non linear equations
%    rJ:        Jacobian
%    mult:      Lagrangian multipliers
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2012 Dynare Team
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


steady_state = [];
params = M.params;
check = 0;
options_.steadystate.nocheck = 1;
nl_func = @(x) dyn_ramsey_static_1(x,M,options_,oo);

if options_.steadystate_flag
    k_inst = [];
    instruments = options_.instruments;
    inst_nbr = size(options_.instruments,1);
    for i = 1:inst_nbr
        k_inst = [k_inst; strmatch(options_.instruments(i,:), ...
                                   M.endo_names,'exact')];
    end
    ys = oo.steady_state;
    if inst_nbr == 1
        inst_val = csolve(nl_func,oo.steady_state(k_inst),'',options_.solve_tolf,100);
    else
        [inst_val,info1] = dynare_solve(nl_func,ys(k_inst),0);
    end
    ys(k_inst) = inst_val;
    exo_ss = [oo.exo_steady_state oo.exo_det_steady_state];
    [xx,params,check] = evaluate_steady_state_file(ys,exo_ss,M,options_);
    [junk,jun,steady_state] = nl_func(inst_val);
else
    n_var = M.orig_endo_nbr;
    xx = oo.steady_state(1:n_var);
    [xx,info1] = dynare_solve(nl_func,xx,0);
    [junk,junk,steady_state] = nl_func(xx);
end



function [resids,rJ,steady_state] = dyn_ramsey_static_1(x,M,options_,oo)
resids = [];
rJ = [];
mult = [];

% recovering usefull fields
params = M.params;
endo_nbr = M.endo_nbr;
endo_names = M.endo_names;
exo_nbr = M.exo_nbr;
orig_endo_nbr = M.orig_endo_nbr;
aux_vars_type = [M.aux_vars.type];
orig_endo_aux_nbr = orig_endo_nbr + min(find(aux_vars_type == 6)) - 1; 
orig_eq_nbr = M.orig_eq_nbr;
inst_nbr = orig_endo_aux_nbr - orig_eq_nbr;
% indices of Lagrange multipliers
i_mult = [orig_endo_aux_nbr+(1:orig_eq_nbr)]';
fname = M.fname;
max_lead = M.maximum_lead;
max_lag = M.maximum_lag;

% indices of all endogenous variables
i_endo = [1:endo_nbr]';
% indices of endogenous variable except instruments
% i_inst = M.instruments;
% lead_lag incidence matrix
i_lag = M.lead_lag_incidence;

if options_.steadystate_flag
    k_inst = [];
    instruments = options_.instruments;
    for i = 1:size(instruments,1)
        k_inst = [k_inst; strmatch(instruments(i,:), ...
                                   endo_names,'exact')];
    end
    oo.steady_state(k_inst) = x;
    [x,params,check] = evaluate_steady_state_file(oo.steady_state,...
                                                  [oo.exo_steady_state; ...
                                                  oo.exo_det_steady_state], ...
                                                  M,options_);
end

xx = zeros(endo_nbr,1);
xx(1:length(x)) = x;
% setting steady state of auxiliary variables
% that depends on original endogenous variables
if any([M.aux_vars.type] ~= 6)
    needs_set_auxiliary_variables = 1;
    fh = str2func([M.fname '_set_auxiliary_variables']);
    s_a_v_func = @(z) fh(z,... 
                         [oo.exo_steady_state,...
                        oo.exo_det_steady_state],...
                         params);
    xx = s_a_v_func(xx);
else
    needs_set_auxiliary_variables = 0;
end

% value and Jacobian of objective function
ex = zeros(1,M.exo_nbr);
[U,Uy,Uyy] = feval([fname '_objective_static'],x,ex, params);
Uy = Uy';
Uyy = reshape(Uyy,endo_nbr,endo_nbr);

% set multipliers and auxiliary variables that
% depends on multipliers to 0 to compute residuals
if (options_.bytecode)
   [chck, res, junk] = bytecode('static',xx,[oo.exo_simul oo.exo_det_simul], ...
               M.params, 'evaluate'); 
   fJ = junk.g1;
else
   [res,fJ] = feval([fname '_static'],xx,[oo.exo_simul oo.exo_det_simul], ...
               M.params);
end
% index of multipliers and corresponding equations
% the auxiliary variables before the Lagrange multipliers are treated
% as ordinary endogenous variables
aux_eq = [1:orig_endo_aux_nbr, orig_endo_aux_nbr+orig_eq_nbr+1:size(fJ,1)];
A = fJ(aux_eq,orig_endo_aux_nbr+1:end);
y = res(aux_eq);
aux_vars = -A\y;

resids1 = y+A*aux_vars;
if inst_nbr == 1
    r1 = sqrt(resids1'*resids1);
else
    [q,r,e] = qr([A y]');
    k = size(A,1)+(1-inst_nbr:0);
    r1 = r(end,k)';
end
if options_.steadystate_flag
    resids = r1;
else
    resids = [res(orig_endo_nbr+(1:orig_endo_nbr-inst_nbr)); r1];
end
rJ = [];
steady_state = [xx(1:orig_endo_aux_nbr); aux_vars];