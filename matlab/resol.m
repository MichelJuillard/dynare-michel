function [dr,info,M,options,oo] = resol(check_flag,M,options,oo)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info},@var{M},@var{options},@var{oo}] =} resol (@var{check_flag},@var{M},@var{options},@var{oo})
%! @anchor{resol}
%! @sp 1
%! Computes first and second order reduced form of the DSGE model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item check_flag
%! Integer scalar, equal to 0 if all the approximation is required, positive if only the eigenvalues are to be computed.
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @item info==30
%! Ergodic variance can't be computed.
%! @end table
%! @sp 1
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{dynare_estimation_init}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! None.
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2011 Dynare Team
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

global it_

jacobian_flag = 0;

if isfield(oo,'dr');
    dr = oo.dr;
end

options = set_default_option(options,'jacobian_flag',1);
info = 0;

it_ = M.maximum_lag + 1 ;

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end

params0 = M.params;

% check if steady_state_0 (-> oo.steady_state) is steady state
tempex = oo.exo_simul;
oo.exo_simul = repmat(oo.exo_steady_state',M.maximum_lag+M.maximum_lead+1,1);
if M.exo_det_nbr > 0
    tempexdet = oo.exo_det_simul;
    oo.exo_det_simul = repmat(oo.exo_det_steady_state',M.maximum_lag+M.maximum_lead+1,1);
end
steady_state = oo.steady_state;
check1 = 0;
% testing for steadystate file
if (~options.bytecode)
    fh = [M.fname '_static'];
end

if options.steadystate_flag
    [steady_state,check1] = feval([M.fname '_steadystate'],steady_state,...
                           [oo.exo_steady_state; ...
                        oo.exo_det_steady_state]);
    if size(steady_state,1) < M.endo_nbr
        if length(M.aux_vars) > 0
            steady_state = add_auxiliary_variables_to_steadystate(steady_state,M.aux_vars,...
                                                           M.fname,...
                                                           oo.exo_steady_state,...
                                                           oo.exo_det_steady_state,...
                                                           M.params,...
                                                           options.bytecode);
        else
            error([M.fname '_steadystate.m doesn''t match the model']);
        end
    end

else
    % testing if steady_state_0  (-> oo.steady_state) isn't a steady state or if we aren't computing Ramsey policy
    if  options.ramsey_policy == 0
        if options.linear == 0
            % nonlinear models
            if (options.block == 0 && options.bytecode == 0)
                if max(abs(feval(fh,steady_state,[oo.exo_steady_state; ...
                                        oo.exo_det_steady_state], M.params))) > options.dynatol
                    [steady_state,check1] = dynare_solve(fh,steady_state,options.jacobian_flag,...
                                                  [oo.exo_steady_state; ...
                                        oo.exo_det_steady_state], M.params);
                end
            else
                [steady_state,check1] = dynare_solve_block_or_bytecode(steady_state,...
                                                                [oo.exo_steady_state; ...
                                    oo.exo_det_steady_state], M.params);
            end;
        else
            if (options.block == 0 && options.bytecode == 0)
                % linear models
                [fvec,jacob] = feval(fh,steady_state,[oo.exo_steady_state;...
                                    oo.exo_det_steady_state], M.params);
                if max(abs(fvec)) > 1e-12
                    steady_state = steady_state-jacob\fvec;
                end
            else
                [steady_state,check1] = dynare_solve_block_or_bytecode(steady_state,...
                                                                [oo.exo_steady_state; ...
                                    oo.exo_det_steady_state], M.params);
            end;
        end
    end
end

% test if M.params_has changed.
if options.steadystate_flag
    updated_params_flag = max(abs(M.params-params0))>1e-12;
else
    updated_params_flag = 0;
end

% testing for problem.
dr.ys = steady_state;

if check1
    if options.steadystate_flag
        info(1)= 19;
        resid = check1 ;
    else
        info(1)= 20;
        resid = feval(fh,oo.steady_state,oo.exo_steady_state, M.params);
    end
    info(2) = resid'*resid ;
    return
end

if ~isreal(steady_state)
    info(1) = 21;
    info(2) = sum(imag(steady_state).^2);
    steady_state = real(steady_state);
    dr.ys = steady_state;
    return
end

if ~isempty(find(isnan(steady_state)))
    info(1) = 22;
    info(2) = NaN;
    dr.ys = steady_state;
    return
end

if options.steadystate_flag && updated_params_flag && ~isreal(M.params)
    info(1) = 23;
    info(2) = sum(imag(M.params).^2);
    dr.ys = steady_state;
    return
end

if options.steadystate_flag && updated_params_flag  && ~isempty(find(isnan(M.params)))
    info(1) = 24;
    info(2) = NaN;
    dr.ys = steady_state;
    return
end


if options.block
    [dr,info,M,options,oo] = dr_block(dr,check_flag,M,options,oo);
else
    [dr,info,M,options,oo] = dr1(dr,check_flag,M,options,oo);
end
if info(1)
    return
end

if M.exo_det_nbr > 0
    oo.exo_det_simul = tempexdet;
end
oo.exo_simul = tempex;
tempex = [];