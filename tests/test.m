function test(a,b,n)
  if max(max(abs(a)-abs(b))) > 1e-5
    error(['Test error in test ' int2str(n)])
  end