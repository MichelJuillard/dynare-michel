function [dr_]=resol(ys,algo,%linear,iorder)
// Copyright (C) 2001 Michel Juillard
// 
 
xlen = xkmax_+xkmin_+1;
klen = ykmin_+ykmax_+1;
iyv = iy_';
iyv = iyv(:);
iyr0 = find(iyv)';
it_ = ykmin_+1;
 
 
if exo_nbr==0 then
  exe_ = [];
end
 
if ~(iy_(ykmin_+1,:)>0) then
  error('RESOL: Error in model specification: some variables don""t appear as current');
end
 
if ykmax_==0 then
  error('Backward or static model: no point in using RESOL');
end
 
if xlen>1 then
  error('RESOL: stochastic exogenous variables must appear only at the'+' current period. Use additional endogenous variables');
end
 
dr_ = struct();
// check if ys is steady state
fff_name = fname_+'_fff';
if max(abs(evstr(fff_name+'(ys)'))) > dynatol_ then
  [ys_temp,%check]=solve(fff_name,ys);
  dr_('ys') = ys_temp;
  if %check then
    error('RESOL: convergence problem in SOLVE');
  end
else
  dr_('ys') = ys
end
 
dr_ = dr1(iorder,dr_);
 
// 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
//            in dr.ghs2 
