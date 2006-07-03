function [e_1,e_2,e_inf]=brm(v,nbr,drop)
e_1=[];e_2=[];e_inf=[];
// Copyright (C) 2001 Michel Juillard
// 
global('y_','iy_','iter_','ykmin_','ykmax_','it_','endo_nbr','ex_');
 
i = (iy_>0)';
iyr0 = find(i(:))';
 
ex_old = ex_;
 
//! mtlb(end) can be replaced by end() or end whether end is an m-file or not
ex_ = ex_(drop+1:mtlb(end),:);
 
//! mtlb(end) can be replaced by end() or end whether end is an m-file or not
y = y_(:,drop+1:mtlb(end));
y = y(:);
 
iter = iter_-drop-ykmin_-ykmax_;
 
z = zeros(endo_nbr,iter);
 
for it_ = ykmin_+1:iter+ykmin_
  z(:,it_) = ff_(y(iyr0));
  iyr0 = iyr0+endo_nbr;
end
 
t1 = z(nbr,:)';
 
//! mtlb(end) can be replaced by end() or end whether end is an m-file or not
t2 = v(drop+1:mtlb(end)-2);
t = t1 ./ t2;
 
//! mtlb_mean(abs(t)) may be replaced by 
//!   mean(abs(t)) if abs(t)is a vector
//!   mean(abs(t),1) if abs(t)is a matrix
e_1 = log(mtlb_mean(abs(t)))/log(10);
 
//!! Unknown function var ,the original calling sequence is used
e_2 = log(var(t))/log(10);
 
//!  mtlb_max(abs(t)) may be replaced by 
//!     max(abs(t)) if abs(t)is a vector
//!     max(abs(t),'r') if abs(t)is a matrix
e_inf = log(mtlb_max(abs(t)))/log(10);
 
ex_ = ex_old;
