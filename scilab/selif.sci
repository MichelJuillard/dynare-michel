function [x]=selif(a,b)
x=[];
// Copyright (C) 2001 Michel Juillard
// 
 
if size(b,2)~=1 then
  error('The second argument in SELIF must be à column vector');
end
 
x = a(matrix(find(b==1),1,-1),:);
 
return
 
 
