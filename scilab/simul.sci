function []=simul(dr,algo,%linear,order,replic)
// Copyright (C) 2001 Michel Juillard
// 
 
global ys_ y_ scalv_ 
 
if algo==[]|algo==0 then
  if ~valf_ then
    if (size(ys_,1)==1)&(ys_==0) then
      ys_ = zeros(size(ys_,1),1);
    end
    y_ = ys_*ones(1,iter_+ykmin_+ykmax_);
    if endval_==1 then
      y_(:,1:ykmin_) = ys0_*ones(1,ykmin_);
    end
  end
   
  if scalv_==[]|scalv_==0 then
    scalv_ = ys_;
  end
   
  scalv_ = 1;
   
  if (size(iy_,2)-nnz(iy_(ykmin_+1,:)))>0 then
    mess = 'DYNARE: error in model specification : variable '+lgy_(matrix(find(iy_(ykmin_+1,:)==0),1,-1),:);
    mess = mess+' doesn''t appear as current variable.';
    error(mess);
  end
   
  if (ykmin_==1)&(ykmax_==1) then
    sim1();
  else
    simk();
  end
else
  simult(dr,replic,iorder,istoch);
end

dyn2vec();

 
// 06/18/01 MJ added dyn2vec if 40 variables or less
// 01/19/03 MJ dyn2vec all endogenous variables whatever their number

 
 
 
 
 
 
