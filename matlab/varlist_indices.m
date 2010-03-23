function [i_var,nvar] = varlist_indices(sublist,list)
% function [i_var,nvar] = varlist_indices(sublist,list)
% returns the indices of a list of endogenous variables
%
% INPUT
%   sublist:    sublist of variables
%   list:       list of variables 
%
% OUTPUT
%   i_var:      variable indices in M_.endo_names
%   nvar:       number of variables in varlist
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2010 Dynare Team
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


[check,i_var] = ismember(sublist,list,'rows');

nvar = length(i_var);

if ~all(check)
    k =find(check);
    error(strcat(sublist(k,:),' hasn''t been declared'))
end

