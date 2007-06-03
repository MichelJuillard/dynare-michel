function writedata(fname)
% function writedata(fname)
% store endogenous and exogenous variables in a text file 
% INPUT
%   fname: name of the text file
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright (2007)
% Gnu Public License.
  global M_ oo_
  S=[fname '_endo.dat'];
  fid = fopen(S,'w');
  for i = 1:size(M_.endo_names,1)
    fprintf(fid,'%s ',M_.endo_names(i,:)');
  end;
  fprintf(fid,'\n');
  for i = 1:size(oo_.endo_simul,2)
    fprintf(fid,'%15.7f ',oo_.endo_simul(:,i));
    fprintf(fid,'\n');
  end
  fclose(fid);
  
  S=[fname '_exo.dat'];
  fid = fopen(S,'w');
  for i = 1:size(M_.exo_names,1)
    fprintf(fid,'%s ',M_.exo_names(i,:));
  end;
  fprintf(fid,'\n');
  for i = 1:size(oo_.exo_simul,1)
    fprintf(fid,'%15.7f ',oo_.exo_simul(i,:));
    fprintf(fid,'\n');
  end
  fclose(fid);
return;