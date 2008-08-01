
function dynasave(s,var_list)

% function dynasave(s,var_list)
% This optional command saves the simulation results in a .MAT file.
% DYNASAVE :	DYNASAVE ( [ 'filename' ] )	
%
% INPUTS
%    s:              matrix of endogenous variables values
%    var_list:       vector of selected endogenous variables
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2007 Dynare Team
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

  global M_ oo_

  n = size(var_list,1);
  if n == 0
    n = M_.endo_nbr;
    ivar = [1:n]';
    var_list = M_.endo_names;
  else
    ivar=zeros(n,1);
    for i=1:n
      i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
	error (['One of the specified variables does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end


%  dyn2vec(var_list(1),var_list(1));
eval([var_list(1) '=oo_.endo_simul(ivar(1),:)'';'])
eval(['save ' s ' ' var_list(1) ' -mat'])
for dynare__i_ = 2:n
    eval([var_list(dynare__i_) '=oo_.endo_simul(ivar(dynare__i_),:)'';'])
    eval(['save ' s ' ' var_list(dynare__i_) ' -append -mat'])
end


