function []=rplot(s1,rplottype)
[nargout,nargin] = argn(0)
// RPLOT :  RPLOT ( ['var1'; 'var2';  ...] , rplottype ) 
// Copyright (C) 2001 Michel Juillard
// 
//  This optionnal command creates the plot of the variable
//  trajectory. By default, the entire simulation period is 
//  ploted. The instruction DSAMPLE permits to reduce the 
//  number of periods in the plot.
 
global dsmpl_ y_

if nargin==1 then
  rplottype = 0;
end
 
ix = (1-ykmin_:size(y_,2)-ykmin_)';
 
m = size(s1,1);
k=[];
for j=1:m
  k1 = grep_exact(lgy_,s1(j));
  if k1==[] then
    error('One of the variable specified does not exist');
  end
  k = [k ; k1];
end
 
 
y = y_(k,:)';
if dsmpl_==0 then
  i = (ykmin_:size(y_,2))';
else
  i = (dsmpl_(1):dsmpl_(2))';
end
 
t = 'Plot of ';
if rplottype==0 then
  xset('window',max(winsid()+1))
  xset('use color',0)
  leg = "";
  for j=1:m
    t = t+s1(j)+' ';
    if j > 1
      leg = leg+'@'
    end
    leg = leg+s1(j);
  end
  plot2d(ix(i),y(i,:),leg=leg);
  xtitle(t,'Periods',' ');
elseif rplottype==1 then
  for j = 1:size(y,1)
    xset('window',max(winsid()+1))
    max(winsid());
    plot2d(ix(i),y(j,i));
    xtitle('Plot of '+s1(:,j),' ',' ');
    xtitle(' ','Periods',' ');
  end
elseif rplottype==2 then
  xset('window',max(winsid()+1))
  max(winsid());
  nl = max(1,fix(size(y,1)/4));
  nc = ceil(size(y,1)/nl);
  for j = 1:size(y,1)
    %v2$2 = int((j-1)/nc);%v2$1 = j-1-nc*%v2$2
    xsetech([%v2$1/nc,%v2$2/nl,1/nc,1/nl]);
    mtlb_plot(ix(i),y(j,i));
    mtlb_hold('on');
    mtlb_plot(ix(i),ys_(j)*ones(1,size(i,1)),'w:');
    xtitle(' ','Periods',' ');
    xtitle(' ',' ',s1(:,j));
    xtitle('Plot of '+s1(:,j),' ',' ');
  end
end
 
// 02/28/01 MJ replaced bseastr by MATLAB's strmatch
// 06/19/01 MJ added 'exact' to strmatch calls
// 09/25/01 MJ translated to Scilab 
 
 
 
 
 
 
 
 
