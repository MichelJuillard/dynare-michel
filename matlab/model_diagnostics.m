function model_diagnostics(M_,options_,oo_)
% function model_diagnostics(M_,options_,oo_)
%   computes various diagnostics on the model 
% INPUTS
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   none
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
    
    endo_names = M_.endo_names;
    lead_lag_incidence = M_.lead_lag_incidence;
    maximum_lag = M_.maximum_lag;
    maximum_lead = M_.maximum_lead;
    
% missing variables at the current period
    k = find(lead_lag_incidence(maximum_lag+1,:)==0);
    if ~isempty(k)
        disp(['The following endogenous variables aren''t present at ' ...
                 'the current period in the model:'])
        for i=1:length(k)
            disp(endo_names(k(i),:))
        end
    end
    