function [z,zss]=dyn2vec(s1,s2)
z=[];zss=[];
[nargout,nargin] = argn(0)
// Copyright (C) 2001 Michel Juillard
// 
 
if dsmpl_==0 then
  %kk = 1:size(y_,2);
else
  %kk = ykmin_+dsmpl_(1):ykmin_+dsmpl_(2);
end
 
if nargin==0 then
  if nargout>1 then
    t = 'DYNARE dyn2vec error: the function doesn''t return values when'+' used without input argument';
    error(t);
  end
  for %i = 1:size(y_,1)
    tt = ['global '+string(lgy_(%i,:)) ; string(lgy_(%i,:))+' = y_(%i,%kk);'];
    deff('[]=assignin()',tt);
    assignin();
    clear('assignin')
  end
  return
else
  j = grep_exact(lgy_,s1);
  if ~(j==[]) then
    %z = y_(j,%kk)';
  else
    j = grep_exact(lgx_,s1);
    if ~(j==[]) then
      if dsmpl_==0 then
        %z = ex_(:,j);
      else
        %z = ex_(xkmin_+dsmpl_(1):xkmin_+dsmpl_(2));
      end
    else
        t = 'DYNARE dyn2vec error: variable '+string(s1(%i,:))+' doesn''t'+' exist.';
      error(t);
    end
  end
end
 
if nargout==1 then
  deff('[]=assignin()',['global '+string(s1) ; string(s1)+' = %z;']);
  assignin();
  clear('assignin')
else
  zss = ys_(j);
end
 
// 02/23/01 MJ redone, incorporating FC's improvements
// 08/24/01 MJ replaced globalize by internal assignin
// 08/24/01 MJ added 'exact' to strmatch (thanks to David Vavra)
// 09/24/01 MJ translated to SciLab 
 
 
 
