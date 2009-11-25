function [i_var,nvar] = varlist_indices(varlist)
% function [i_var,nvar] = varlist_indices(varlist)
% returns the indices of a list of endogenous variables
%
% INPUT
%   varlist:    (character area) list of variables
%
% OUTPUT
%   i_var:      variable indices in M_.endo_names
%   nvar:       number of variables in varlist
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009 Dynare Team
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

    global M_
    
    endo_nbr = M_.endo_nbr;
    
    if isempty(varlist)
        varlist = M_.endo_names(1:M_.orig_endo_nbr,:);
    end
        i_var = [];
        for i=1:size(varlist,1)
            tmp = strmatch(varlist(i,:),M_.endo_names,'exact');
            if isempty(tmp)
                error([tmp ' isn''t an endogenous variable'])
            end
            i_var = [i_var; tmp];
        end
        nvar = length(i_var);
