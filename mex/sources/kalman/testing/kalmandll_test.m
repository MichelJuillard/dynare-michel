function [LIKDLL loglik]=kalmandll_test(T,mf,R,Q,H,Pstar,Pinf,data,start)

if isempty(H)
    H=zeros(size(data,1), size(data,1))
elseif H==0
    H=zeros(size(data,1), size(data,1))
end
Z=zeros(size(data,1), size(T,2))
for i = 1:size(data,1)
Z(i,mf(i))=1
end
LIKDLL= kalman_filter_dll4(T,Z,R,Q,H,Pstar,data,start)
%Y=data;
if isempty(Pinf)
    Pinf=zeros(size(T));
elseif Pinf==0
    Pinf=zeros(size(T));
end
% test DiffuseLikelihoodH1
loglik = dynare_filter(Z,H,T,R,Q,data,Pstar,Pinf)