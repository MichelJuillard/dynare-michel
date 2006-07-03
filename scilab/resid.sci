function []=resid()
// Copyright (C) 2001 Michel Juillard
// 
global('iter_','valf_','ex_','y_','it_','exe_','ys_','iy_','ykmin_','ykmax_','endval_','z');
 
n = size(iy_,2);
//  if ~ valf_ | size(y_,2) ~= iter_+ykmin_+ykmax_
if ~valf_ then
  if (size(ys_,1)==1)&(ys_==0) then
    ys_ = zeros(size(ys_,1),1);
  end
  y_ = ys_*ones(1,iter_+ykmin_+ykmax_);
  if endval_==1 then
     
    //! mtlb(ys0_) can be replaced by ys0_() or ys0_ whether ys0_ is an m-file or not
    y_(:,1:ykmin_) = mtlb(ys0_)*ones(1,ykmin_);
  end
end
 
i = iy_';
iyr0 = find(i(:))';
 
y = y_(:);
z = zeros(n,iter_);
for it_ = ykmin_+1:iter_+ykmin_
  z(:,it_-ykmin_) = ff_(y(iyr0));
  iyr0 = iyr0+n;
end
 
dyn_disp([(1:iter_)',z']);
 
 
