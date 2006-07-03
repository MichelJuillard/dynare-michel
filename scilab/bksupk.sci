function [d1]=bksupk(ny,fid,jcf,icc1)
d1=[];
// Copyright (C) 2001 Michel Juillard
// 
 
global('ykmax_','iter_','c','ncc');
 
icf = 1:jcf-1;
ir = (iter_-1)*ny+1:ny*iter_;
%irf = icc1+(iter_-1)*ny;
d1 = zeros(iter_*ny,1);
 
ofs = ((iter_-1)*ny+1-1)*ncc*8;
junk = mseek(ofs,fid,'set');
c = mtlb_fread(fid,[ncc,ny],'float64');
c = c';
 
%v = c(:,jcf)
d1(ir,1) = %v(:);
ir = ir-ny;
 
i = 2;
 
while i<=ykmax_|i<=iter_ then
  irf1 = selif(%irf,%irf<iter_*ny);
   
  ofs = ((iter_-i)*ny+1-1)*ncc*8;
  junk = mseek(ofs,fid,'set');
  c = mtlb_fread(fid,[ncc,ny],'float64');
  c = c';
   
  %v1 = c(:,jcf)-c(:,1:size(irf1,1))*d1(irf1)
  d1(ir,1) = %v1(:);
  ir = ir-ny;
  %irf = %irf-ny;
  i = i+1;
end
 
while i<=iter_ then
   
  ofs = ((iter_-i)*ny+1-1)*ncc*8;
  junk = mseek(ofs,fid,'set');
  c = mtlb_fread(fid,[ncc,ny],'float64');
  c = c';
   
  %v1 = c(:,jcf)-c(:,icf)*d1(%irf)
  d1(ir,1) = %v1(:);
  ir = ir-ny;
  %irf = %irf-ny;
  i = i+1;
end

 
