function []=linear(x)
// Copyright (C) 2001 Michel Juillard
// 
// LINEAR : LINEAR (x,'filename')
 
 
x = matrix(x,max(size(x)),1);
xx = x;
nn = size(iy_,2);
nv = max(iy_(size(iy_,1),:));
indic = zeros(0,0);
 
if size(x,1)==nv then
  r = matrix(find(matrix(iy_'>0,(ykmin_+ykmax_+1)*nn,1)),1,-1);
  n = abs(lgy_);
  n = ones(ykmin_+ykmax_+1,1).*.n;
  n = ascii(n(r,:));
  m = ((-ykmin_:ykmax_)')*ones(1,nn);
  m = matrix(m',(ykmin_+ykmax_+1)*nn,1);
  m = m(r);
elseif size(x,1)==nn then
  l = iy_>0;
  l = bool2s(l) .* (ones(size(iy_,1),1).*.(1:nn));
   
  //!! Unknown function nonzeros ,the original calling sequence is used
  i = nonzeros(matrix(iy_,size(iy_,1)*nn,1));
   
  //!! Unknown function nonzeros ,the original calling sequence is used
  j = nonzeros(matrix(l,size(iy_,1)*nn,1));
  s = ones(1,max(iy_(size(iy_,1),:)));
  indic = sparse([i(:),j(:)],s,[nv,nn]);
  x = indic*x;
  n = lgy_;
  m = zeros(nn,1);
else
  error('Wrong number of arguments in LINEAR.');
end
 
jacob(x,'ff_');
 
if ~(indic==[]) then
  jacobia_ = jacobia_*indic;
  clear('indic');
end
 
mtlb_fprintf(1,'Periods  :  ');
mtlb_fprintf(1,'%4g \n',iter_);
mtlb_fprintf(1,'Endogenous variables : ');
mtlb_fprintf(1,'\n');
mtlb_fprintf(1,lgy_);
mtlb_fprintf(1,'\n');
mtlb_fprintf(1,'Exogenous variables : ');
mtlb_fprintf(1,'\n');
mtlb_fprintf(1,lgx_);
mtlb_fprintf(1,'\n');
mtlb_fprintf(1,'Linearization around :');
mtlb_fprintf(1,'\n');
 
for i = 1:size(n,1)
  mtlb_fprintf(1,n(i,:));
  mtlb_fprintf(1,'(%1g)',m(i));
  mtlb_fprintf(1,' = %15.6f \n',xx(i));
end
 
mtlb_fprintf(1,'\n');
 
for i = 1:size(jacobia_,1)
  for j = 1:size(jacobia_,2)
    if jacobia_(i,j)~=0 then
      if jacobia_(i,j)==1 then
        if j==1 then
          mtlb_fprintf(1,n(j,:));
          mtlb_fprintf(1,'(%1g)',m(j));
        else
          mtlb_fprintf(1,' + ');
          mtlb_fprintf(1,n(j,:));
          mtlb_fprintf(1,'(%1g)',m(j));
        end
      elseif jacobia_(i,j)==(-1) then
        mtlb_fprintf(1,' - ');
        mtlb_fprintf(1,n(j,:));
        mtlb_fprintf(1,'(%1g)',m(j));
      elseif jacobia_(i,j)>0 then
        if j>1 then
          mtlb_fprintf(1,' + ');
        end
        mtlb_fprintf(1,'%15.6g',jacobia_(i,j));
        mtlb_fprintf(1,'*');
        mtlb_fprintf(1,n(j,:));
        mtlb_fprintf(1,'(%1g)',m(j));
      else
        mtlb_fprintf(1,'%15.6g',jacobia_(i,j));
        mtlb_fprintf(1,'*');
        mtlb_fprintf(1,n(j,:));
        mtlb_fprintf(1,'(%1g)',m(j));
         
      end
    end
  end
  mtlb_fprintf(1,'\n');
end
 
return
 
