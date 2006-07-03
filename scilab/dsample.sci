function []=dsample(s1,s2)
[nargout,nargin] = argn(0)
// Copyright (C) 2001 Michel Juillard
// 
// DSAMPLE : DSAMPLE(d1,d2)
//  This optional command permits to reduce the number of
//  periods considered in following output commands. 
//               If only one argument is 
//  provided, output is from period -maxlag+1 to the period 
//  specified in the DSAMPLE command. If two arguments are
//  present output is done for the interval between the 
//  two periods.
//               DSAMPLE without arguments reset the sample to the one 
//               specified by PERIODS 
 
global dsmpl_
 
dsmpl_ = zeros(2,1);
 
if nargin==0 then
  dsmpl_(1) = 1
  dsmpl_(2,1) = iter_(:)+ykmin_+ykmax_;
elseif nargin==1 then
  dsmpl_(1) = 1
  dsmpl_(2,1) = s1(:)+ykmin_;
else
  dsmpl_(1,1) = s1(:)+ykmin_;
  dsmpl_(2,1) = s2(:)+ykmin_;
end
 
if dsmpl_(1) > iter_+ykmin_|dsmpl_(2) > iter_+ykmin_+ykmax_ then
  t = 'DYNARE dsample error: one of the arguments is larger than the one'+' specified in PERIODS';
  error(t);
end
 
// 02/23/01 MJ added error checking
// 02/19/01 MJ added ykmin_


