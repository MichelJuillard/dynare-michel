function [x] = shift(x, varargin)
%% mode 1 : [ x1(t), x1(t-1), ..., x2(t), x2(t-1), ... ]
%% mode 2 : [ x1(t), x2(t), ..., x1(t-1), x2(t-1), ... ]

z = -1;
mode = 2;

if ( nargin > 1 )
  z = varargin(1)';
end

if ( nargin > 2 )
  mode = varargin(2);
end

[n, k] = size(x);
nz = length(z);

if isempty(z)
  x = zeros([n, 0]);
  return
end

maxz = max([0, z]);
minz = abs(min([0, z]));
x = [ nan*zeros([minz,k]); x; nan*zeros([maxz,k]) ];
%n = n + maxz + minz;

m = (1 : n)';
m = m(:, ones([1, k*nz]));

switch mode
case 1
  t = 0 : (n+maxz+minz) : (k-1)*(n+maxz+minz);
  t = t(ones([nz,1]),:);
  t = t(:)';
  t = t(ones([n,1]),:);
  z = z(ones([n,1]),:);
  z = z(:);
  z = reshape(z(:, ones([1,k])), n, k*nz);
case 2
  t = ( 0 : (n+maxz+minz) : (k-1)*(n+maxz+minz) )';
  t = t(:,ones([nz,1]));
  t = t(:)';
  t = t(ones([n,1]),:);
  z = z(ones([k,1]),:);
  z = z(:)';
  z = z(ones([n,1]),:);
end  

m = m + t + z;
m = m + minz;
x = x(:);
x = x(m);

return
