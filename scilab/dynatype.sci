function []=dynatype(s)
// Copyright (C) 2001 Michel Juillard
// 
// DYNATYPE : DYNATYPE ( [ 'filename' ] )
//  This optional command saves the simulation
//  results in a text file. The name of each 
//  variable preceeds the corresponding results.
//  This command must follow SIMUL.
 
[fid,%v] = mopen(s,'w',0)
if %v<0 then fid = -1;end
for i = 1:size(lgy_,1)
  mtlb_fprintf(fid,lgy_(i,:),'\n');
  mtlb_fprintf(fid,'\n');
  mtlb_fprintf(fid,'%15.8g\n',y_(i,:)');
end
mclose(fid);
 
return
 
