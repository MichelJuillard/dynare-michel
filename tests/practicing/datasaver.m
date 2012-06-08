function datasaver (s,var_list)
% datasaver saves variables simulated by Dynare
% INPUT
%    s: a string containing the name of the destination *.m file
%    var_list: a character matrix containting the name of the variables
%              to be saved (optional, default: all endogenous variables)
% OUTPUT
%    none
% This is part of the examples included in F. Barillas, R. Colacito,
% S. Kitao, C. Matthes, T. Sargent and Y. Shin (2007) "Practicing
% Dynare".
  
% Modified by M. Juillard to make it also compatible with Dynare 
% version 4 (12/4/07)   
  

global lgy_ lgx_ y_ endo_nbr M_ oo_

% test and adapt for Dynare version 4
if isempty(lgy_)
  lgy_ = M_.endo_names;
  lgx_ + M_.exo_names;
  y_ = oo_.endo_simul;
  endo_nbr = M_.endo_nbr;
end

sm=[s,'.m'];
fid=fopen(sm,'w') ;

n = size(var_list,1);
if n == 0
  n = endo_nbr;
  ivar = [1:n]';
  var_list = lgy_;
else
  ivar=zeros(n,1);
  for i=1:n
    i_tmp = strmatch(var_list(i,:),lgy_,'exact');
    if isempty(i_tmp)
      error (['One of the specified variables does not exist']) ;
    else
      ivar(i) = i_tmp;
    end
  end
end


for i = 1:n
	fprintf(fid,[lgy_(ivar(i),:), '=['],'\n') ;
	fprintf(fid,'\n') ;
	fprintf(fid,'%15.8g\n',y_(ivar(i),:)') ;
   	fprintf(fid,'\n') ;
   	fprintf(fid,'];\n') ;
  	fprintf(fid,'\n') ;
end
fclose(fid) ;

return ;
