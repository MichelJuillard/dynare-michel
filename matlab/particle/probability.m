function [prior,likelihood,C,posterior] = probability(mu,sqrtP,prior,X)
[dim,nov] = size(X);              
M = size(mu,2) ;
if nargout>1
  likelihood = zeros(M,nov);        
  normfact = (2*pi)^(dim/2);  
  for k=1:M
    XX = bsxfun(@minus,X,mu(:,k));
    S = sqrtP(:,:,k);
    foo = S \ XX;
    likelihood(k,:) = exp(-0.5*sum(foo.*foo, 1))/abs((normfact*prod(diag(S))));
  end
end
likelihood = likelihood + 1e-99;
if nargout>2
  C = prior*likelihood + 1e-99;                   
end
if nargout>3
  posterior = bsxfun(@rdivide,bsxfun(@times,prior',likelihood),C) + 1e-99 ;
  posterior = bsxfun(@rdivide,posterior,sum(posterior,1));
end
