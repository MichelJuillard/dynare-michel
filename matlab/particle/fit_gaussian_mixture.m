function [StateMu,StateSqrtP,StateWeights] = fit_gaussian_mixture(X,StateMu,StateSqrtP,StateWeights,crit,niters,check) 
[dim,Ndata] = size(X);             
M = size(StateMu,2) ;
if check                        % Ensure that covariances don't collapse
  MIN_COVAR_SQRT = sqrt(eps);
  init_covars = StateSqrtP;
end
eold = -Inf;
for n=1:niters
  % Calculate posteriors based on old parameters
  [prior,likelihood,marginal,posterior] = probability(StateMu,StateSqrtP,StateWeights,X);
  e = sum(log(marginal));
  if (n > 1 && abs((e - eold)/eold) < crit)
    return;
  else
    eold = e;
  end
  new_pr = (sum(posterior,2))';
  StateWeights = new_pr/Ndata;
  StateMu = bsxfun(@rdivide,(posterior*X')',new_pr);
  for j=1:M
    diffs = bsxfun(@minus,X,StateMu(:,j));
    tpost = (1/sqrt(new_pr(j)))*sqrt(posterior(j,:));
    diffs = bsxfun(@times,diffs,tpost);
    [foo,tcov] = qr2(diffs',0);
    StateSqrtP(:,:,j) = tcov';
    if check
      if min(abs(diag(StateSqrtP(:,:,j)))) < MIN_COVAR_SQRT
        StateSqrtP(:,:,j) = init_covars(:,:,j);
      end
    end
  end
end     

