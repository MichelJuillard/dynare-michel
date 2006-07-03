function [y]=reshapel(x,m,n)
y=[];
// Copyright (C) 2001 Michel Juillard
// 
//RESHAPEL Change size.
// RESHAPEL(X,M,N) returns the M-by-N matrix whose elements
// are taken linewise from X.  An error results if X does
// not have M*N elements.
 
[mm,nn] = size(x);
 
if mm*nn~=m*n then
  error('Matrix must have M*N elements.');
end
 
%v = x';
[%i,%i] = find(%v);%i=%i(:);%i=%i(:);
%v = %v(:)
if %i<>[] then s = %v(%i+size(%v,1)*(%i-1)) ;else s = [],end
if size(i,2)~=1 then
  i = i';
end
k = (j-1)*nn+i;
j = k-1-fix(k-1./n).*n+1;
i = (k-j)/n+1;
y = full(sparse([i(:),j(:)],s,[m,n]));
 
return
 
