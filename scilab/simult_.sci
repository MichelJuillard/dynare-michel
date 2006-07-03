function [y_]=simult_(y0,dr,ex_,iorder,iter)
y_=[];
// Copyright (C) 2001 Michel Juillard
// 
// 

y_ = zeros(endo_nbr,iter+ykmin_);
y_(:,1:ykmin_) = ones(1,ykmin_).*.y0;
if iorder==1 then
  for i = ykmin_+1:iter+ykmin_
    temp = y_(dr('order_var'),i-1)-ones(1,ykmin_).*.dr('ys')(dr('order_var'));
    temp = temp(dr('nstatic')+1:dr('nstatic')+dr('npred'));
    y_(dr('order_var'),i) = dr('ys')(dr('order_var'))+dr('ghx')*temp+dr('ghu')*ex_(i+xkmin_-ykmin_,:)';
  end
elseif iorder==2 then
  for i = ykmin_+1:iter+ykmin_
    tempx = y_(dr('order_var'),i-1)-ones(1,ykmin_).*.dr('ys')(dr('order_var'));
    tempx = tempx(dr('nstatic')+1:dr('nstatic')+dr('npred'));
    tempu = ex_(i+xkmin_-ykmin_,:)';
    tempxx = tempx.*.tempx;
    tempxu = tempx.*.tempu;
    tempuu = tempu.*.tempu;
    y_(dr('order_var'),i) = dr('ys')(dr('order_var'))+dr('ghs2')/2+dr('ghx')*tempx+dr('ghu')*tempu+.5*(dr('ghxx')*tempxx+dr('ghuu')*tempuu)+dr('ghxu')*tempxu;
  end
end
 
