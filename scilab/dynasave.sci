function []=dynasave(s)
// Copyright (C) 2001 Michel Juillard
// 
// DYNASAVE : DYNASAVE ( [ 'filename' ] ) 
//  This optional command saves the simulation results
//  in a .BIN datafile. This command must follow SIMUL.
//  Filename must be specified without the .BIN extension.
 
 
if ~(matrix(find(abs(s)==46),1,-1)==[]) then
  dyn_disp('Warning : Error in DYNASAVE');
  dyn_disp('Filename '+s+' is incorrect.');
  dyn_disp('You must enter filename without extension,');
  dyn_disp('For example : '+s(1,1:matrix(find(abs(s)==46),1,-1)-1));
  return
   
end
 
s = s+'.BIN';
 
[fid,%v] = mopen(s,'w',0)
if %v<0 then fid = -1;end
mtlb_fwrite(fid,size(y_),'int');
mtlb_fwrite(fid,size(lgy_,2),'int');
mtlb_fwrite(fid,ykmin_,'int');
mtlb_fwrite(fid,ykmax_,'int');
mtlb_fwrite(fid,xkmin_,'int');
mtlb_fwrite(fid,xkmax_,'int');
mtlb_fwrite(fid,lgy_,'int');
mtlb_fwrite(fid,y_,'float64');
mclose(fid);
 
return
 


