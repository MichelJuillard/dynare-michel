function [a]=indnv(x,y)
a=[];
// Copyright (C) 2001 Michel Juillard
// 
 
%v=size(x)
a = zeros(%v(1),%v(2));
 
for i = 1:size(x,1)
  j = find(x(i)==y)';
  if j==[] then
    a(i) = %nan;
  else
    a(i) = j(1);
  end
end
 
 
 
