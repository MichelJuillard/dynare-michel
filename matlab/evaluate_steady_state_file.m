function [ys,params1,check] = evaluate_steady_state_file(ys_init,exo_ss,params,fname,steadystate_flag)
% function [ys,params1,check] = evaluate_steady_state_file(ys_init,exo_ss,params,steadystate_flag)
% Evaluates steady state files 
%  
% INPUTS
%   ys_init                   vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state
%   params                    vector           parameters
%   steadystate_check_flag    boolean          if true, check that the
%                                              steadystate verifies the
%                                              static model         
%  
% OUTPUTS
%   ys                        vector           steady state
%   params1                   vector           model parameters possibly
%                                              modified by user steadystate
%                                              function
%   check                     2x1 vector       error codes
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


    if steadystate_flag == 1
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
    