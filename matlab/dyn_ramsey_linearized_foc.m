function [jacobia_,dr,info,M_,oo_] = dyn_ramsey_linearized_foc(dr,M_,options_,oo_)
% function [jacobia_,dr,info,M_,oo_] = dyn_ramsey_linearized_foc(dr,M_,options_,oo_)
% computes the Jacobian of the linear approximation of the F.O.C of a
% Ramsey problem
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code. 
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%                                 info=6: The jacobian matrix evaluated at the steady state is complex.        
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

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

info = 0;

if isfield(M_,'orig_model')
        orig_model = M_.orig_model;
        M_.endo_nbr = orig_model.endo_nbr;
        M_.orig_endo_nbr = orig_model.orig_endo_nbr;
        M_.aux_vars = orig_model.aux_vars;
        M_.endo_names = orig_model.endo_names;
        M_.lead_lag_incidence = orig_model.lead_lag_incidence;
        M_.maximum_lead = orig_model.maximum_lead;
        M_.maximum_endo_lead = orig_model.maximum_endo_lead;
        M_.maximum_lag = orig_model.maximum_lag;
        M_.maximum_endo_lag = orig_model.maximum_endo_lag;
    end

    if options_.steadystate_flag
        k_inst = [];
        instruments = options_.instruments;
        for i = 1:size(instruments,1)
            k_inst = [k_inst; strmatch(options_.instruments(i,:), ...
                                       M_.endo_names,'exact')];
        end
        ys = oo_.steady_state;
        [inst_val,info1] = dynare_solve('dyn_ramsey_static_', ...
                                oo_.steady_state(k_inst),0, ...
                                M_,options_,oo_,it_);
        M_.params = evalin('base','M_.params;');
        ys(k_inst) = inst_val;
        [x,check] = feval([M_.fname '_steadystate'],...
                          ys,[oo_.exo_steady_state; ...
                            oo_.exo_det_steady_state]);
        if size(x,1) < M_.endo_nbr 
            if length(M_.aux_vars) > 0
                x = add_auxiliary_variables_to_steadystate(x,M_.aux_vars,...
                                                           M_.fname,...
                                                           oo_.exo_steady_state,...
                                                           oo_.exo_det_steady_state,...
                                                           M_.params);
            else
                error([M_.fname '_steadystate.m doesn''t match the model']);
            end
        end
        oo_.steady_state = x;
        [junk,junk,multbar] = dyn_ramsey_static_(oo_.steady_state(k_inst),M_,options_,oo_,it_);
    else
        [oo_.steady_state,info1] = dynare_solve('dyn_ramsey_static_', ...
                                        oo_.steady_state,0,M_,options_,oo_,it_);
        [junk,junk,multbar] = dyn_ramsey_static_(oo_.steady_state,M_,options_,oo_,it_);
    end
        
    check1 = max(abs(feval([M_.fname '_static'],...
                           oo_.steady_state,...
                           [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state], M_.params))) > options_.dynatol ;
    if check1
        info(1) = 20;
        info(2) = check1'*check1;
        return
    end
    
    [jacobia_,M_] = dyn_ramsey_dynamic_(oo_.steady_state,multbar,M_,options_,oo_,it_);
    klen = M_.maximum_lag + M_.maximum_lead + 1;
    dr.ys = [oo_.steady_state;zeros(M_.exo_nbr,1);multbar];
    oo_.steady_state = dr.ys;
    
    if options_.noprint == 0
        disp_steady_state(M_,oo_)
        for i=M_.orig_endo_nbr:M_.endo_nbr
            if strmatch('mult_',M_.endo_names(i,:))
                disp(sprintf('%s \t\t %g',M_.endo_names(i,:), ...
                             dr.ys(i)));
            end
        end
    end

