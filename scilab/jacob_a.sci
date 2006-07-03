function [jacobia_]=jacob_a(func,z,varargin)
// Copyright (C) 2001 Michel Juillard
// 
// symmetric formula to compute the Jacobian
// 
 
global d1_
nargin = length(varargin);

if nargin > 2 then
  d1_ = evstr(func+'(z,varargin)');
else
  d1_ = evstr(func+'(z)');
end
nz = size(z,1);
jacobia_ = zeros(size(d1_,1),nz);
dh = max(abs(z),gstep_*ones(nz,1))*(%eps^(1/3));
xdh1 = z;
xdh2 = z;
for j = 1:nz
  xdh1(j) = z(j)-dh(j);
  xdh2(j) = z(j)+dh(j);
  h = xdh2(j)-xdh1(j);
  if nargin > 2 then
    jacobia_(:,j) = (evstr(func+'(xdh2,varargin)')-evstr(func+'(xdh1,varargin)'))/h;
  else
    jacobia_(:,j) = (evstr(func+'(xdh2)')-evstr(func+'(xdh1)'))/h;
  end
  xdh1(j) = z(j);
  xdh2(j) = z(j);
end
 
// 10/21/02 MJ creation
