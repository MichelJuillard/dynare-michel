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
%LIKDLL= kalman_filter_dll8(T,Z,R,Q,H,Pstar,data,start)
%[loglik,per,d] = kalman_filter_dll(Z,H,T,R,Q,data,a,Pstar)
[LIK2 lik2] = kalman_filter(T,R,Q,H,Pstar,data,start,mf,options_.kalman_tol,options_.riccati_tol)
if isempty(Pinf)
    Pinf=zeros(size(T));
elseif Pinf==0
    Pinf=zeros(size(T));
end
% test DiffuseLikelihoodH1
%[loglikd1,per,d] = kalman_filter_dll(Z,H,T,R,Q,data,a,Pstar,Pinf)
[LIKdlikd ] = diffuse_kalman_filter(T,R,Q,H,Pinf,Pstar,data,start,Z,options_.kalman_tol,options_.riccati_tol)
%loglikd2 = dynare_filter(Z,H,T,R,Q,data,Pstar,Pinf)
