function [y]=fx_(x)
y=[];
// Copyright (C) 2001 Michel Juillard
// 
global('ys_','ykmin_','ykmax_','it_','iy_','xkmin_','ex_');
 
//!! Unknown function repmat ,the original calling sequence is used
y = repmat(ys_,ykmin_+ykmax_+1,1);
iyr0 = matrix(find(matrix(iy_',size(iy_,1)*size(iy_,2),1)),1,-1);
y = y(iyr0);
ex_(it_+xkmin_-ykmin_,:) = x';
y = ff_(y);
