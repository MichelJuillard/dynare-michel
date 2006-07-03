function [x]=shift(x,varargin)
[nargout,nargin] = argn(0)
//% mode 1 : [ x1(t), x1(t-1), ..., x2(t), x2(t-1), ... ]
//% mode 2 : [ x1(t), x2(t), ..., x1(t-1), x2(t-1), ... ]
 
z = -1;
%mode = 2;
 
if nargin>1 then
  z = varargin(1)';
end
 
if nargin>2 then
  %mode = varargin(2);
end
 
[n,k] = size(x);
 
//! unknown arg type, using mtlb_length 
nz = mtlb_length(z);
 
if z==[] then
  %v1=[n,0]
  x = zeros(%v1(1),%v1(2));
  return
   
end
 
maxz = max([0,z]);
minz = abs(min([0,z]));
%v=[minz,k]
%v=[maxz,k]
x = [%nan*zeros(%v(1),%v(2));x;%nan*zeros(%v(1),%v(2))];
//n = n + maxz + minz;
 
m = (1:n)';
%v=[1,k*nz]
m = m(:,ones(%v(1),%v(2)));
 
select %mode
case 1 then
  t = 0:n+maxz+minz:(k-1)*(n+maxz+minz);
  %v1=[nz,1]
  t = t(ones(%v1(1),%v1(2)),:);
  t = t(:)';
  %v1=[n,1]
  t = t(ones(%v1(1),%v1(2)),:);
  %v1=[n,1]
  z = z(ones(%v1(1),%v1(2)),:);
  z = z(:);
  %v1=[1,k]
  z = matrix(z(:,ones(%v1(1),%v1(2))),n,k*nz);
case 2 then
  t = (0:n+maxz+minz:(k-1)*(n+maxz+minz))';
  %v1=[nz,1]
  t = t(:,ones(%v1(1),%v1(2)));
  t = t(:)';
  %v1=[n,1]
  t = t(ones(%v1(1),%v1(2)),:);
  %v1=[k,1]
  z = z(ones(%v1(1),%v1(2)),:);
  z = z(:)';
  %v1=[n,1]
  z = z(ones(%v1(1),%v1(2)),:);
end
 
m = m+t+z;
m = m+minz;
x = x(:);
x = x(m);
 
return
 
