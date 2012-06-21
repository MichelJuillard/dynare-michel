function [ys,params,info] = evaluate_steady_state_file(ys_init,exo_ss,M,options)
% function [ys,params1,info] = evaluate_steady_state_file(ys_init,exo_ss,M,options)
% Evaluates steady state files 
%  
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state
%   M                         struct           model parameters
%   options                   struct           options
%  
% OUTPUTS
%   ys                        vector           steady state
%   params1                   vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   info                      2x1 vector       error codes
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2012 Dynare Team
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

    ys = [];
    params = [];
    info = 0; 
    params = M.params;
    fname = M.fname;
    if options.steadystate_flag == 1
        % old format
        assignin('base','tmp_00_',params);
        evalin('base','M_.params=tmp_00_; clear(''tmp_00_'')');
        h_steadystate = str2func([fname '_steadystate']);                       
        [ys,check] = h_steadystate(ys_init, exo_ss);
        params1 = evalin('base','M_.params');
    else % steadystate_flag == 2
         % new format
        h_steadystate = str2func([fname '_steadystate2']);                       
        [ys,params1,check] = h_steadystate(ys_init, exo_ss, params);
    end            
    
    if check
        info(1) = 19
        info(2) = NaN;
    end
    
    updated_params_flag = max(abs(params1-params)) > 1e-12 || ~isequal(isnan(params1),isnan(params)); %checks whether numbers or NaN changed

    h_set_auxiliary_variables = str2func([M.fname '_set_auxiliary_variables']);
    if  isnan(updated_params_flag) || (updated_params_flag  && any(isnan(params(~isnan(params))-params1(~isnan(params))))) %checks if new NaNs were added
        info(1) = 24;
        info(2) = NaN;
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
        return
    end

    if updated_params_flag && ~isreal(params1)
        info(1) = 23;
        info(2) = sum(imag(params).^2);
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
        return
    end
    
    if updated_params_flag
        params = params1;
    end

    % adding values for auxiliary variables
    if length(M.aux_vars) > 0 && ~options.ramsey_policy
        ys = h_set_auxiliary_variables(ys,exo_ss,params);
    end

    check1 = 0;
    if ~options.steadystate.nocheck
        % Check whether the steady state obtained from the _steadystate file is a steady state.
        [residuals,check] = evaluate_static_model(ys,exo_ss,params,M,options);
        if check
            info(1) = 19;
            info(2) = check; % to be improved
            return;
        end
        if max(abs(residuals)) > options.dynatol.f
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

    
