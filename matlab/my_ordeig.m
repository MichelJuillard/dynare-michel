function eval = my_ordeig(t)
  
  n = size(t,2);
  eval = zeros(n,1);
  for i=1:n-1
    if t(i+1,i) == 0
      eval(i) = t(i,i);
    else
      k = i:i+1;
      eval(k) = eig(t(k,k));
      i = i+1;
    end
  end
  if i < n
    t(n) = t(n,n);
  end
  
      