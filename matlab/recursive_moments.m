function [mu,sigma] = recursive_moments(m0,s0,data,offset)
% stephane.adjemian@cepremap.cnrs.fr [02-03-2005]
%
% m0    :: the initial mean, assumed to be a n*1 vector.
% s0    :: the initial variance, assumed to be a n*n matrix.
% data  :: T*n matrix.  
%
% mu    :: the sample mean, a n*1 vector.
% sigma :: the sample covariance matrix, a n*n matrix.
%
% --> TODO: dll  
T = size(data,1);
n = size(data,2);

for t = 1:T
  tt = t+offset;
  m1 = m0 + (1/tt)*(data(t,:)'-m0);
  s1 = s0 + (1/tt)*(data(t,:)'*data(t,:)-m1*m1'-s0) + ((tt-1)/tt)*(m0*m0'-m1*m1');
  m0 = m1;
  s0 = s1;
end

mu = m1;
sigma = s1;


  

