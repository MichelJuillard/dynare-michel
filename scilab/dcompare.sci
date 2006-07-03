function []=dcompare(s1)
// Copyright (C) 2001 Michel Juillard
// 
 
global('dsmpl_','iter_');
global('nvx','nvy','x','y','lag1');
 
ftest(s1,0);
 
i = (lag1(1):size(x,2)-lag1(2)+1)';
 
if size(dsmpl_,1)==1 then
  error('DSAMPLE not specified.');
end
 
if dsmpl_(3)>0 then
  if dsmpl_(3)==2 then
    if dsmpl_(1)<0|dsmpl_(2)>(size(x,2)-lag1(2)) then
      error('Wrong sample.');
    end
    i = (dsmpl_(1)+lag1(1):dsmpl_(2)+lag1(1))';
  elseif dsmpl_(3)==1 then
    if dsmpl_(1)>(size(x,2)-lag1(2)) then
      error('Wrong sample.');
    end
    i = (lag1(1):dsmpl_(1)+lag1(1))';
  end
end
 
j = bseastr(nvx,nvy);
 
 
//! mtlb(stop) can be replaced by stop() or stop whether stop is an m-file or not
if mtlb(stop) then
  return
   
end
 
 
//! mtlb_mean(abs(x(j,i)-y(j,i))) may be replaced by 
//!   mean(abs(x(j,i)-y(j,i))) if abs(x(j,i)-y(j,i))is a vector
//!   mean(abs(x(j,i)-y(j,i)),1) if abs(x(j,i)-y(j,i))is a matrix
z = mean(mtlb_mean(abs(x(j,i)-y(j,i))));
 
dyn_disp('The mean absolute difference between set '+s1(1,:)+'and set '+s1(2,:));
dyn_disp('is : '+string(z));
return
 
 
 
