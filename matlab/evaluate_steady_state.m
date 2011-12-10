function [ys,params,info] = evaluate_steady_state(ys_init,M,options,oo,steadystate_check_flag)
% function [ys,info] = evaluate_steady_state(M,options,oo)
% Computes the steady state 
%  
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   M                         struct           model structure
%   options                   struct           options
%   oo                        struct           output results
%   steadystate_check_flag    boolean          if true, check that the
%                                              steadystate verifies the
%                                              static model         
%  
% OUTPUTS
%   ys                        vector           steady state
%   params                    vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

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

    info = 0;
    check = 0;
    
    steadystate_flag = options.steadystate_flag;
    params = M.params;
    exo_ss = [oo.exo_steady_state; oo.exo_det_steady_state];
    updated_params_flag = 0;
    
    if length(M.aux_vars) > 0
        h_set_auxiliary_variables = str2func([M.fname '_set_auxiliary_variables']);                       
        ys_init = h_set_auxiliary_variables(ys_init,exo_ss,M.params);
    end
    
    if options.ramsey_policy
        [ys,params] = dyn_ramsey_static(ys_init,M,options,oo);
    elseif steadystate_flag
        % explicit steady state file
        [ys,params1,check] = evaluate_steady_state_file(ys_init,exo_ss,params,M.fname,steadystate_flag);
        updated_params_flag = max(abs(params1-params)) > 1e-12;
        if updated_params_flag
            params = params1;
        end
        % adding values for auxiliary variables
        if length(M.aux_vars) > 0
            ys = h_set_auxiliary_variables(ys,exo_ss,M.params);
        end
        check1 = 0;
        if steadystate_check_flag
            % Check whether the steady state obtained from the _steadystate file is a steady state.
            [residuals,check] = evaluate_static_model(ys,exo_ss,params,M,options);
            if check
                info(1) = 19;
                info(2) = check; % to be improved
                return;
            end
            if max(abs(residuals)) > options.dynatol
            info(1) = 19;
            info(2) = residuals'*residuals;
            return
        end
    elseif ~isempty(options.steadystate_partial)
            ssvar = options.steadystate_partial.ssvar;
            nov   = length(ssvar);
            indv  = zeros(nov,1);
            for i = 1:nov
                indv(i) = strmatch(ssvar(i),M.endo_names,'exact');
            end
            [ys,check] = dynare_solve('restricted_steadystate',...
                                                    ys(indv),...
                                                    options.jacobian_flag, ...         
                                                    exo_ss,indv);
        end
    elseif (options.bytecode == 0 && options.block == 0)
        if options.linear == 0
            % non linear model
            [ys,check] = dynare_solve([M.fname '_static'],...
                                      ys_init,...
                                      options.jacobian_flag, ...     
                                      exo_ss, params);
        else
            % linear model
            fh_static = str2func([M.fname '_static']);
            [fvec,jacob] = fh_static(ys_init,exo_ss, ...
                                 params);
            if max(abs(fvec)) > 1e-12
                ys = ys_init-jacob\fvec;
            else
                ys = ys_init;
            end

        end
    else
        % block or bytecode
        [ys,check] = dynare_solve_block_or_bytecode(ys_init,exo_ss, params, options, M);
    end

    if check
        if options.steadystate_flag
            info(1)= 19;
            resid = check1 ;
        else
            info(1)= 20;
            resid = evaluate_static_model(ys_init,exo_ss,params,M,options);
        end
        info(2) = resid'*resid ;
        return
    end

    if ~isreal(ys)
        info(1) = 21;
        info(2) = sum(imag(ys).^2);
        ys = real(ys);
        return
    end

    if ~isempty(find(isnan(ys)))
        info(1) = 22;
        info(2) = NaN;
        return
    end

    if options.steadystate_flag && updated_params_flag && ~isreal(M.params)
        info(1) = 23;
        info(2) = sum(imag(M.params).^2);
        return
    end

    if options.steadystate_flag && updated_params_flag  && ~isempty(find(isnan(M.params)))
        info(1) = 24;
        info(2) = NaN;
        return
    end

