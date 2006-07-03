function [x]=fff(y)
x=[];
// Copyright (C) 2001 Michel Juillard
// 
 
global('ykmin_','ykmax_','iy_');
 
iyr = matrix(find(matrix(iy_'>0,(ykmin_+ykmax_+1)*size(iy_,2),1)),1,-1);
y = ones(ykmin_+ykmax_+1,1).*.y;
y = y(iyr,:);
x = ff_(y);
 
return
 
