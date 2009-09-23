function [JJ,H] = mc_getJ(lpmat, M_,oo_,options_,0,indx,indexo,mf2,nlags,useautocorr);

for j=1:size(lpmat,1),
set_parameters(lpmat(j,:));
[JJ, H] = getJJ(M_,oo_,options_,0,indx,indexo,mf2,nlags,useautocorr);
dynare_identification(H,JJ);
% [indok, indno, indweak] = dynare_identification(H,JJ);
end
