
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
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.


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


