function []=steady()
// Copyright (C) 2001 Michel Juillard
// 
 
steady_();
for i = 1:size(ys_,1)
  mprintf('%s \t\t %g\n',lgy_(i),ys_(i));
  mfprintf(fh_log,'%s \t\t %g\n',lgy_(i),ys_(i));
end
 
// 06/24/01 MJ steady print results; steady_ doesn't
// 02/06/02 MJ chanded to dyn_mprintf()
