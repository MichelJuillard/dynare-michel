% Copyright (C) 2001 Michel Juillard
%
function dynatype (s,var_list)
% DYNATYPE :	DYNATYPE ( [ 'filename' ] )
%		This optional command saves the simulation
%		results in a text file. The name of each 
%		variable preceeds the corresponding results.
%		This command must follow SIMUL.

global M_ oo_

fid=fopen(s,'w') ;

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


for i = 1:n
	fprintf(fid,M_.endo_names(ivar(i),:),'\n') ;
	fprintf(fid,'\n') ;
	fprintf(fid,'%15.8g\n',oo_.y_simul(ivar(i),:)') ;
end
fclose(fid) ;

return ;
