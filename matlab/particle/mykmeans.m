function [c,SqrtVariance,Weights] = mykmeans(x,g,init,cod) 
[n,m] = size(x) ;
indold = zeros(1,m) ;
if cod==0
  d = transpose(sum(bsxfun(@power,bsxfun(@minus,x,mean(x)),2)));
  d = sortrows( [transpose(1:m) d],2) ;
  d = d((1+(0:1:g-1))*m/g,1) ;
  c = x(:,d);
else
  c = init ;
end 
for iter=1:300 
  dist = zeros(g,m) ;
  for i=1:g
    dist(i,:) = sum(bsxfun(@power,bsxfun(@minus,x,c(:,i)),2));
  end
  [rien,ind] = min(dist) ;
  if isequal(ind,indold) 
    break ;
  end
  indold = ind ;
  for i=1:g 
    lin = bsxfun(@eq,ind,i.*ones(1,m)) ;
    h = x(:,lin) ;
    c(:,i) = mean(h,2) ;
  end
end
SqrtVariance = zeros(n,n,g) ; 
Weights = zeros(1,g) ; 
for i=1:g
  temp = x(:,bsxfun(@eq,ind,i*ones(1,m))) ;
  u = bsxfun(@minus,temp,mean(temp,2)); %temp-mean(temp,1)' ;
  SqrtVariance(:,:,i) = chol( (u*u')/size(temp,2) )' ;
  Weights(i) = size(temp,2)/m ;
end